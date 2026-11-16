//! Independent Rust replay validator for source bundles, variant collections, and development bundles.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::Read;
use std::path::Path;

#[cfg(unix)]
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};

use crate::crypto::{base64_decode_strict, sha256_hex, sha512t24u};
use crate::gzip::decompress_gzip;
use crate::jcs::{
    canonical_bytes, digest_object_without, digest_value, jcs_string, parse_json_object,
    parse_json_object_with_limit, JsonValue,
};

pub const MAX_JSON_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_SOURCE_DESCRIPTOR_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_GENBANK_ARTIFACT_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_SEQUENCE_INPUT_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_GFF3_JSON_MEMBERS: usize = 5_000_000;
pub const MAX_VARIANTS: usize = 100_000;
pub const MAX_VARIANT_BASES: usize = 1_000_000;
pub const MAX_IDENTIFIER_BYTES: usize = 256;

const RAW_PROFILE: &str = "raw-iupac-dna/v1";
const FASTA_PROFILE: &str = "fasta-dna/v1";
const REFERENCE_FASTA_PROFILE: &str = "fasta-reference-dna/v1";
const GENBANK_PROFILE: &str = "genbank-273-traditional-dna-physical-structural/v2";
const GFF3_PROFILE: &str = "gff3-external-sequence/v1";
const EXTERNAL_PROFILE: &str = "external-dna-source/v1";

#[cfg(target_os = "macos")]
const OPEN_NOFOLLOW_NONBLOCK: i32 = 0x0100 | 0x0004;
#[cfg(all(target_os = "linux", any(target_arch = "x86_64", target_arch = "x86")))]
const OPEN_NOFOLLOW_NONBLOCK: i32 = 0o400000 | 0o4000;
#[cfg(all(target_os = "linux", any(target_arch = "aarch64", target_arch = "arm")))]
const OPEN_NOFOLLOW_NONBLOCK: i32 = 0o100000 | 0o4000;
#[cfg(all(unix, not(any(
    target_os = "macos",
    all(target_os = "linux", any(target_arch = "x86_64", target_arch = "x86", target_arch = "aarch64", target_arch = "arm"))
))))]
const OPEN_NOFOLLOW_NONBLOCK: i32 = 0;

fn is_sha256_hex(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

fn is_b64u_32(s: &str) -> bool {
    s.len() == 32
        && s.bytes()
            .all(|b| matches!(b, b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_'))
}

fn is_refget_id(s: &str) -> bool {
    s.starts_with("SQ.") && is_b64u_32(&s[3..])
}

fn is_iupac_base(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

pub fn read_regular_file(path: &Path, label: &str, max_bytes: usize) -> Result<Vec<u8>, String> {
    let link_meta =
        fs::symlink_metadata(path).map_err(|e| format!("cannot stat {label}: {e}"))?;
    if link_meta.file_type().is_symlink() || !link_meta.is_file() {
        return Err(format!("{label} must be a non-symlink regular file"));
    }
    #[cfg(unix)]
    if link_meta.nlink() != 1 {
        return Err(format!(
            "{label} must have link count 1 (found {})",
            link_meta.nlink()
        ));
    }
    if link_meta.len() > max_bytes as u64 {
        return Err(format!("{label} exceeds maximum byte size {max_bytes}"));
    }
    let mut opts = OpenOptions::new();
    opts.read(true);
    #[cfg(unix)]
    if OPEN_NOFOLLOW_NONBLOCK != 0 {
        opts.custom_flags(OPEN_NOFOLLOW_NONBLOCK);
    }
    let mut file: File = opts.open(path).map_err(|e| format!("cannot open {label}: {e}"))?;
    let fd_meta = file
        .metadata()
        .map_err(|e| format!("cannot fstat {label}: {e}"))?;
    if !fd_meta.is_file() {
        return Err(format!("{label} is not a regular file"));
    }
    #[cfg(unix)]
    {
        if fd_meta.nlink() != 1 {
            return Err(format!("{label} must have link count 1"));
        }
        if (fd_meta.dev(), fd_meta.ino()) != (link_meta.dev(), link_meta.ino()) {
            return Err(format!("{label} changed identity during open"));
        }
    }
    let mut buf = Vec::new();
    let mut limited = file.by_ref().take((max_bytes + 1) as u64);
    limited
        .read_to_end(&mut buf)
        .map_err(|e| format!("cannot read {label}: {e}"))?;
    if buf.len() > max_bytes {
        return Err(format!("{label} exceeded maximum byte size during read"));
    }
    let after_meta = file
        .metadata()
        .map_err(|e| format!("cannot fstat {label} after read: {e}"))?;
    if buf.len() as u64 != fd_meta.len() || after_meta.len() != fd_meta.len() {
        return Err(format!("{label} size changed during read"));
    }
    #[cfg(unix)]
    if (after_meta.dev(), after_meta.ino(), after_meta.nlink(), after_meta.mtime_nsec(), after_meta.ctime_nsec())
        != (fd_meta.dev(), fd_meta.ino(), fd_meta.nlink(), fd_meta.mtime_nsec(), fd_meta.ctime_nsec())
    {
        return Err(format!("{label} changed during read"));
    }
    Ok(buf)
}

fn read_dir_files(
    dir: &Path,
    role_limits: &BTreeMap<&str, usize>,
    required: &[&str],
    label: &str,
) -> Result<BTreeMap<String, Vec<u8>>, String> {
    let meta = fs::symlink_metadata(dir).map_err(|e| format!("cannot stat {label}: {e}"))?;
    if meta.file_type().is_symlink() || !meta.is_dir() {
        return Err(format!("{label} must be a non-symlink directory"));
    }
    let mut out = BTreeMap::new();
    for entry_res in fs::read_dir(dir).map_err(|e| format!("cannot list {label}: {e}"))? {
        let entry = entry_res.map_err(|e| format!("cannot read {label} entry: {e}"))?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| format!("non-UTF8 filename in {label}"))?;
        let Some(&limit) = role_limits.get(name.as_str()) else {
            return Err(format!("unexpected child {name:?} in {label}"));
        };
        let child_bytes = read_regular_file(&entry.path(), &format!("{label}/{name}"), limit)?;
        out.insert(name, child_bytes);
    }
    let after = fs::symlink_metadata(dir).map_err(|e| format!("cannot restat {label}: {e}"))?;
    if after.file_type().is_symlink() || !after.is_dir() {
        return Err(format!("{label} changed type during read"));
    }
    #[cfg(unix)]
    if (after.dev(), after.ino()) != (meta.dev(), meta.ino()) {
        return Err(format!("{label} changed identity during read"));
    }
    for req in required {
        if !out.contains_key(*req) {
            return Err(format!("{label} is missing required child {req:?}"));
        }
    }
    Ok(out)
}

fn expect_keys<'a>(
    val: &'a JsonValue,
    expected: &[&str],
    label: &str,
) -> Result<&'a [(String, JsonValue)], String> {
    let entries = val
        .as_object()
        .ok_or_else(|| format!("{label} must be an object"))?;
    let expected_set: BTreeSet<&str> = expected.iter().copied().collect();
    let actual_set: BTreeSet<&str> = entries.iter().map(|(k, _)| k.as_str()).collect();
    if actual_set != expected_set {
        return Err(format!("{label} keys invalid"));
    }
    Ok(entries)
}

fn verify_artifact_seal(obj: &JsonValue, label: &str) -> Result<String, String> {
    let entries = obj
        .as_object()
        .ok_or_else(|| format!("{label} must be an object"))?;
    let stored = obj
        .get("artifact_sha256")
        .and_then(|v| v.as_str())
        .ok_or_else(|| format!("{label} missing artifact_sha256"))?;
    if !is_sha256_hex(stored) {
        return Err(format!("{label}.artifact_sha256 must be lowercase hex"));
    }
    let computed = digest_object_without(entries, "artifact_sha256")?;
    if computed != stored {
        return Err(format!("{label} artifact_sha256 mismatch"));
    }
    Ok(stored.to_string())
}

fn decode_fasta_bytes(raw: &[u8], wrapper: &str) -> Result<BTreeMap<String, String>, String> {
    let logical: Vec<u8> = match wrapper {
        "gzip" => {
            if !raw.starts_with(b"\x1f\x8b") {
                return Err("gzip FASTA input missing gzip header".to_string());
            }
            decompress_gzip(raw)?
        }
        "identity" => {
            if raw.starts_with(b"\x1f\x8b") {
                return Err("identity FASTA input has unexpected gzip header".to_string());
            }
            raw.to_vec()
        }
        _ => return Err(format!("unsupported FASTA wrapper: {wrapper}")),
    };
    if !logical.is_ascii() {
        return Err("FASTA input must be ASCII".to_string());
    }
    let text = std::str::from_utf8(&logical).map_err(|e| format!("invalid FASTA UTF-8: {e}"))?;
    let mut records: BTreeMap<String, String> = BTreeMap::new();
    let mut current_id: Option<String> = None;
    let mut chunks = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(prev) = current_id.take() {
                records.insert(prev, std::mem::take(&mut chunks));
            }
            let header = rest.trim();
            if header.is_empty() {
                return Err("FASTA record has empty header".to_string());
            }
            let rid = header
                .split_whitespace()
                .next()
                .ok_or_else(|| "FASTA record has empty id".to_string())?
                .to_string();
            if records.contains_key(&rid) {
                return Err(format!("duplicate FASTA record id: {rid}"));
            }
            current_id = Some(rid);
        } else {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if current_id.is_none() {
                return Err("FASTA sequence appears before header".to_string());
            }
            let upper = trimmed.to_ascii_uppercase();
            if !upper.bytes().all(is_iupac_base) {
                return Err("FASTA record contains non-IUPAC bases".to_string());
            }
            chunks.push_str(&upper);
        }
    }
    if let Some(last) = current_id {
        records.insert(last, chunks);
    }
    if records.is_empty() {
        return Err("FASTA input contains no records".to_string());
    }
    Ok(records)
}

fn verify_source_bundle_inner(
    source_dir: &Path,
    source_inputs: &BTreeMap<String, String>,
) -> Result<(JsonValue, BTreeMap<String, JsonValue>, BTreeMap<String, String>), String> {
    let mut role_limits = BTreeMap::new();
    role_limits.insert("source.json", MAX_SOURCE_DESCRIPTOR_BYTES);
    role_limits.insert("sequence.json", MAX_JSON_BYTES);
    role_limits.insert("reference.json", MAX_JSON_BYTES);
    role_limits.insert("external.json", MAX_SOURCE_DESCRIPTOR_BYTES);
    role_limits.insert("genbank.json", MAX_GENBANK_ARTIFACT_BYTES);
    role_limits.insert("annotation.json", MAX_SOURCE_DESCRIPTOR_BYTES);

    let files = read_dir_files(source_dir, &role_limits, &["source.json"], "source bundle")?;
    let desc = parse_json_object(&files["source.json"], "source.json")?;
    expect_keys(
        &desc,
        &[
            "format",
            "version",
            "producer",
            "source_ir",
            "source_ir_sha256",
            "artifact_sha256",
        ],
        "source.json",
    )?;
    if desc.get("format").and_then(|v| v.as_str()) != Some("brainc.source-descriptor")
        || desc.get("version").and_then(|v| v.as_i64()) != Some(1)
    {
        return Err("source.json format/version invalid".to_string());
    }

    let expected_producer = JsonValue::Object(vec![
        ("name".to_string(), JsonValue::String("brainc-source".to_string())),
        ("version".to_string(), JsonValue::String("1.0.0".to_string())),
        (
            "passes".to_string(),
            JsonValue::Array(vec![
                JsonValue::String("select-source-profile".to_string()),
                JsonValue::String("compile-native-source".to_string()),
                JsonValue::String("bind-source-evidence".to_string()),
                JsonValue::String("emit-source-descriptor".to_string()),
            ]),
        ),
    ]);
    let producer = desc
        .get("producer")
        .ok_or_else(|| "source.json missing producer".to_string())?;
    if !producer.exact_eq(&expected_producer) {
        return Err("source.json producer mismatch".to_string());
    }

    let source_ir = desc
        .get("source_ir")
        .ok_or_else(|| "source.json missing source_ir".to_string())?;
    expect_keys(
        source_ir,
        &["profile", "parameters", "inputs", "artifacts", "records"],
        "source_ir",
    )?;
    let stored_ir_sha = desc
        .get("source_ir_sha256")
        .and_then(|v| v.as_str())
        .ok_or_else(|| "source.json missing source_ir_sha256".to_string())?;
    if digest_value(source_ir)? != stored_ir_sha {
        return Err("source_ir_sha256 mismatch".to_string());
    }
    verify_artifact_seal(&desc, "source.json")?;

    let profile = source_ir
        .get("profile")
        .and_then(|v| v.as_str())
        .ok_or_else(|| "source_ir.profile must be a string".to_string())?;
    let expected_children: Vec<(&str, &str, &str, i64, &str)> = match profile {
        RAW_PROFILE | FASTA_PROFILE => vec![(
            "sequence",
            "sequence.json",
            "brain01.sequence-collection-ir",
            1,
            "collection_ir_sha256",
        )],
        REFERENCE_FASTA_PROFILE => vec![(
            "reference",
            "reference.json",
            "brainc.reference-sequence-catalog",
            1,
            "catalog_sha256",
        )],
        GENBANK_PROFILE => vec![(
            "genbank",
            "genbank.json",
            "brainc.bio.insdc-genbank-ir",
            2,
            "bio_ir_sha256",
        )],
        GFF3_PROFILE => vec![
            (
                "sequence",
                "sequence.json",
                "brain01.sequence-collection-ir",
                1,
                "collection_ir_sha256",
            ),
            (
                "annotation",
                "annotation.json",
                "brainc.bio.gff3-ir",
                1,
                "bio_ir_sha256",
            ),
        ],
        EXTERNAL_PROFILE => vec![(
            "external",
            "external.json",
            "brainc.external-source-closure",
            1,
            "closure_ir_sha256",
        )],
        _ => return Err(format!("unsupported source profile: {profile}")),
    };

    let mut expected_file_names: BTreeSet<String> = BTreeSet::new();
    expected_file_names.insert("source.json".to_string());
    for &(_, fname, _, _, _) in &expected_children {
        expected_file_names.insert(fname.to_string());
    }
    let actual_file_names: BTreeSet<String> = files.keys().cloned().collect();
    if actual_file_names != expected_file_names {
        return Err("source bundle filename closure does not match profile".to_string());
    }

    let artifacts_meta = source_ir
        .get("artifacts")
        .and_then(|v| v.as_object())
        .ok_or_else(|| "source_ir.artifacts must be an object".to_string())?;
    if artifacts_meta.len() != expected_children.len() {
        return Err("source_ir.artifacts role count mismatch".to_string());
    }

    let mut parsed_children = BTreeMap::new();
    for &(role, fname, exp_fmt, exp_ver, ir_key) in &expected_children {
        let child_ref = source_ir
            .get("artifacts")
            .and_then(|a| a.get(role))
            .ok_or_else(|| format!("source_ir.artifacts missing {role}"))?;
        expect_keys(
            child_ref,
            &["format", "version", "artifact_sha256", "ir_sha256"],
            &format!("source_ir.artifacts.{role}"),
        )?;
        let limit = if role == "annotation" {
            MAX_GFF3_JSON_MEMBERS
        } else {
            crate::jcs::MAX_JSON_MEMBERS
        };
        let child_obj = parse_json_object_with_limit(&files[fname], fname, limit)?;
        let child_art_sha = verify_artifact_seal(&child_obj, fname)?;
        if child_obj.get("format").and_then(|v| v.as_str()) != Some(exp_fmt)
            || child_obj.get("version").and_then(|v| v.as_i64()) != Some(exp_ver)
        {
            return Err(format!("{fname} format/version mismatch"));
        }
        if child_ref.get("format").and_then(|v| v.as_str()) != Some(exp_fmt)
            || child_ref.get("version").and_then(|v| v.as_i64()) != Some(exp_ver)
            || child_ref.get("artifact_sha256").and_then(|v| v.as_str()) != Some(child_art_sha.as_str())
            || child_ref.get("ir_sha256").and_then(|v| v.as_str())
                != child_obj.get(ir_key).and_then(|v| v.as_str())
        {
            return Err(format!("source_ir.artifacts.{role} disagrees with {fname}"));
        }
        parsed_children.insert(role.to_string(), child_obj);
    }

    let mut sequences = BTreeMap::new();
    match profile {
        RAW_PROFILE | FASTA_PROFILE | GFF3_PROFILE => {
            let seq_art = &parsed_children["sequence"];
            let members = seq_art
                .get("members")
                .and_then(|v| v.as_array())
                .ok_or_else(|| "sequence.json missing members array".to_string())?;
            for m in members {
                let seq_ir = m
                    .get("sequence_artifact")
                    .and_then(|sa| sa.get("sequence_ir"))
                    .ok_or_else(|| "sequence member missing sequence_ir".to_string())?;
                let rid = seq_ir
                    .get("record_id")
                    .and_then(|v| v.as_str())
                    .ok_or_else(|| "sequence_ir missing record_id".to_string())?;
                let seq = seq_ir
                    .get("sequence")
                    .and_then(|v| v.as_str())
                    .ok_or_else(|| "sequence_ir missing sequence".to_string())?;
                let upper = seq.to_ascii_uppercase();
                if sequences.insert(rid.to_string(), upper).is_some() {
                    return Err(format!("duplicate record_id {rid:?} in sequence.json"));
                }
            }
        }
        GENBANK_PROFILE => {
            let gb_art = &parsed_children["genbank"];
            let members = gb_art
                .get("sequence_collection")
                .and_then(|sc| sc.get("members"))
                .and_then(|v| v.as_array())
                .ok_or_else(|| "genbank.json missing sequence_collection.members".to_string())?;
            for m in members {
                let rid = m
                    .get("record_id")
                    .and_then(|v| v.as_str())
                    .ok_or_else(|| "genbank member missing record_id".to_string())?;
                let chunks = m
                    .get("sequence")
                    .and_then(|s| s.get("storage"))
                    .and_then(|st| st.get("chunks"))
                    .and_then(|c| c.as_array())
                    .ok_or_else(|| "genbank member missing sequence chunks".to_string())?;
                let mut seq = String::new();
                for ch in chunks {
                    let part = ch
                        .as_str()
                        .ok_or_else(|| "genbank sequence chunk must be a string".to_string())?;
                    seq.push_str(part);
                }
                let upper = seq.to_ascii_uppercase();
                if sequences.insert(rid.to_string(), upper).is_some() {
                    return Err(format!("duplicate record_id {rid:?} in genbank.json"));
                }
            }
        }
        REFERENCE_FASTA_PROFILE => {
            let seq_path_str = source_inputs.get("sequence").ok_or_else(|| {
                "fasta-reference-dna/v1 requires --source-input sequence=PATH".to_string()
            })?;
            let raw = read_regular_file(
                Path::new(seq_path_str),
                "sequence input",
                MAX_SEQUENCE_INPUT_BYTES,
            )?;
            let input_ref = source_ir
                .get("inputs")
                .and_then(|i| i.get("sequence"))
                .ok_or_else(|| "source_ir.inputs missing sequence".to_string())?;
            if input_ref.get("sha256").and_then(|v| v.as_str()) != Some(sha256_hex(&raw).as_str())
                || input_ref.get("byte_length").and_then(|v| v.as_i64()) != Some(raw.len() as i64)
            {
                return Err("reference sequence input bytes do not match source.json".to_string());
            }
            let wrapper = source_ir
                .get("parameters")
                .and_then(|p| p.get("wrapper"))
                .and_then(|w| w.as_str())
                .ok_or_else(|| "fasta-reference-dna/v1 missing parameters.wrapper".to_string())?;
            sequences = decode_fasta_bytes(&raw, wrapper)?;
        }
        EXTERNAL_PROFILE => {
            let seq_path_str = source_inputs.get("sequence").ok_or_else(|| {
                "external-dna-source/v1 requires --source-input sequence=PATH".to_string()
            })?;
            let raw = read_regular_file(
                Path::new(seq_path_str),
                "sequence input",
                MAX_SEQUENCE_INPUT_BYTES,
            )?;
            sequences = decode_fasta_bytes(&raw, "identity")?;
        }
        _ => {}
    }

    // Verify any explicitly supplied physical source_inputs against source_ir.inputs
    if !matches!(profile, REFERENCE_FASTA_PROFILE | EXTERNAL_PROFILE) && !source_inputs.is_empty() {
        let inputs_map = source_ir
            .get("inputs")
            .and_then(|v| v.as_object())
            .ok_or_else(|| "source_ir.inputs must be an object".to_string())?;
        for (role, path_str) in source_inputs {
            let (_, inp_ref) = inputs_map
                .iter()
                .find(|(k, _)| k == role)
                .ok_or_else(|| format!("unexpected --source-input role {role:?}"))?;
            let raw = read_regular_file(
                Path::new(path_str),
                &format!("source input {role}"),
                MAX_SOURCE_DESCRIPTOR_BYTES,
            )?;
            if inp_ref.get("sha256").and_then(|v| v.as_str()) != Some(sha256_hex(&raw).as_str())
                || inp_ref.get("byte_length").and_then(|v| v.as_i64()) != Some(raw.len() as i64)
            {
                return Err(format!("source input {role:?} bytes disagree with source.json"));
            }
        }
    }

    let records = source_ir
        .get("records")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "source_ir.records must be an array".to_string())?;
    let mut expected_ids = BTreeSet::new();
    for rec in records {
        expect_keys(
            rec,
            &["record_id", "bases", "sequence_sha256", "refget_id"],
            "source_ir.records[]",
        )?;
        let rid = rec
            .get("record_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| "record_id must be string".to_string())?;
        if !expected_ids.insert(rid.to_string()) {
            return Err(format!("duplicate record_id {rid:?} in source_ir.records"));
        }
        let bases = rec
            .get("bases")
            .and_then(|v| v.as_i64())
            .ok_or_else(|| "bases must be int".to_string())?;
        let seq_sha = rec
            .get("sequence_sha256")
            .and_then(|v| v.as_str())
            .ok_or_else(|| "sequence_sha256 must be string".to_string())?;
        let refget_id = rec
            .get("refget_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| "refget_id must be string".to_string())?;
        let seq = sequences
            .get(rid)
            .ok_or_else(|| format!("missing extracted sequence for record {rid:?}"))?;
        if !seq.is_ascii() || !seq.bytes().all(is_iupac_base) {
            return Err(format!("record {rid:?} contains non-IUPAC or non-ASCII bytes"));
        }
        if seq.len() as i64 != bases
            || sha256_hex(seq.as_bytes()) != seq_sha
            || format!("SQ.{}", sha512t24u(seq.as_bytes())) != refget_id
        {
            return Err(format!("record {rid:?} sequence integrity mismatch"));
        }
    }
    let actual_ids: BTreeSet<String> = sequences.keys().cloned().collect();
    if actual_ids != expected_ids {
        return Err("extracted sequence record ids differ from source_ir.records".to_string());
    }

    Ok((desc, parsed_children, sequences))
}

pub fn vrs_location_digest(refget_id: &str, start: i64, end: i64) -> Result<String, String> {
    if !is_refget_id(refget_id) || start < 0 || end < start {
        return Err("invalid VRS SequenceLocation inputs".to_string());
    }
    let obj = JsonValue::Object(vec![
        ("end".to_string(), JsonValue::Int(end)),
        (
            "sequenceReference".to_string(),
            JsonValue::Object(vec![
                (
                    "refgetAccession".to_string(),
                    JsonValue::String(refget_id.to_string()),
                ),
                (
                    "type".to_string(),
                    JsonValue::String("SequenceReference".to_string()),
                ),
            ]),
        ),
        ("start".to_string(), JsonValue::Int(start)),
        (
            "type".to_string(),
            JsonValue::String("SequenceLocation".to_string()),
        ),
    ]);
    Ok(sha512t24u(&canonical_bytes(&obj)?))
}

pub fn vrs_allele_digest(loc_digest: &str, alt: &str) -> Result<String, String> {
    if !is_b64u_32(loc_digest) || !alt.bytes().all(is_iupac_base) {
        return Err("invalid VRS Allele inputs".to_string());
    }
    let obj = JsonValue::Object(vec![
        (
            "location".to_string(),
            JsonValue::String(loc_digest.to_string()),
        ),
        (
            "state".to_string(),
            JsonValue::Object(vec![
                ("sequence".to_string(), JsonValue::String(alt.to_string())),
                (
                    "type".to_string(),
                    JsonValue::String("LiteralSequenceExpression".to_string()),
                ),
            ]),
        ),
        ("type".to_string(), JsonValue::String("Allele".to_string())),
    ]);
    Ok(sha512t24u(&canonical_bytes(&obj)?))
}

pub fn validate_variant_collection_paths(
    source_dir: &Path,
    variant_path: &Path,
    source_inputs: &BTreeMap<String, String>,
) -> Result<JsonValue, String> {
    let mut source_sha_seen: Option<String> = None;
    let mut variant_sha_seen: Option<String> = None;
    let result = validate_variant_collection_inner(
        source_dir,
        variant_path,
        source_inputs,
        &mut source_sha_seen,
        &mut variant_sha_seen,
    );
    match result {
        Ok(report) => Ok(report),
        Err(detail) => {
            let core = JsonValue::Object(vec![
                (
                    "format".to_string(),
                    JsonValue::String("brainc.variant-validation-report".to_string()),
                ),
                ("version".to_string(), JsonValue::Int(1)),
                (
                    "validator".to_string(),
                    JsonValue::Object(vec![
                        (
                            "name".to_string(),
                            JsonValue::String("epigenesis-validate-rs".to_string()),
                        ),
                        ("version".to_string(), JsonValue::String("1.0.0".to_string())),
                    ]),
                ),
                ("valid".to_string(), JsonValue::Bool(false)),
                ("code".to_string(), JsonValue::String("VARVAL001".to_string())),
                (
                    "detail".to_string(),
                    JsonValue::String(format!("VARVAL001: {detail}")),
                ),
                (
                    "source_artifact_sha256".to_string(),
                    source_sha_seen.map(JsonValue::String).unwrap_or(JsonValue::Null),
                ),
                (
                    "variant_artifact_sha256".to_string(),
                    variant_sha_seen.map(JsonValue::String).unwrap_or(JsonValue::Null),
                ),
                (
                    "checks".to_string(),
                    JsonValue::Array(vec![JsonValue::Object(vec![
                        (
                            "name".to_string(),
                            JsonValue::String("variant-replay".to_string()),
                        ),
                        ("passed".to_string(), JsonValue::Bool(false)),
                    ])]),
                ),
                ("summary".to_string(), JsonValue::Null),
            ]);
            let report_sha = digest_value(&core)?;
            let mut entries = core.as_object().unwrap().to_vec();
            entries.push(("report_sha256".to_string(), JsonValue::String(report_sha)));
            Err(String::from_utf8(canonical_bytes(&JsonValue::Object(entries))?).unwrap())
        }
    }
}

fn validate_variant_collection_inner(
    source_dir: &Path,
    variant_path: &Path,
    source_inputs: &BTreeMap<String, String>,
    source_sha_seen: &mut Option<String>,
    variant_sha_seen: &mut Option<String>,
) -> Result<JsonValue, String> {
    let (source_desc, _children, sequences) = verify_source_bundle_inner(source_dir, source_inputs)?;
    let source_art_sha = source_desc
        .get("artifact_sha256")
        .and_then(|v| v.as_str())
        .unwrap()
        .to_string();
    *source_sha_seen = Some(source_art_sha.clone());

    let raw_var = read_regular_file(variant_path, "variant collection", MAX_JSON_BYTES)?;
    let payload = parse_json_object(&raw_var, "variant collection")?;
    if let Some(s) = payload.get("artifact_sha256").and_then(|v| v.as_str()) {
        *variant_sha_seen = Some(s.to_string());
    }
    expect_keys(
        &payload,
        &[
            "format",
            "version",
            "vrs_version",
            "producer",
            "source",
            "variants_ir",
            "variants_ir_sha256",
            "artifact_sha256",
        ],
        "variant collection",
    )?;
    if payload.get("format").and_then(|v| v.as_str()) != Some("brainc.variant-collection")
        || payload.get("version").and_then(|v| v.as_i64()) != Some(1)
        || payload.get("vrs_version").and_then(|v| v.as_str()) != Some("2.x")
    {
        return Err("invalid variant collection format/version/vrs_version".to_string());
    }

    let expected_producer = JsonValue::Object(vec![
        ("name".to_string(), JsonValue::String("brainc-variant".to_string())),
        ("version".to_string(), JsonValue::String("1.0.0".to_string())),
        (
            "passes".to_string(),
            JsonValue::Array(vec![
                JsonValue::String("validate-source-bundle".to_string()),
                JsonValue::String("extract-reference-sequences".to_string()),
                JsonValue::String("normalize-and-verify-variants".to_string()),
                JsonValue::String("compute-vrs-v2-identifiers".to_string()),
                JsonValue::String("seal-variant-collection".to_string()),
            ]),
        ),
    ]);
    if !payload.get("producer").unwrap().exact_eq(&expected_producer) {
        return Err("unexpected variant producer object".to_string());
    }

    let src_meta = payload.get("source").unwrap();
    expect_keys(
        src_meta,
        &[
            "format",
            "version",
            "profile",
            "source_ir_sha256",
            "artifact_sha256",
            "records",
        ],
        "variant.source",
    )?;
    let source_ir = source_desc.get("source_ir").unwrap();
    if src_meta.get("format").and_then(|v| v.as_str()) != source_desc.get("format").and_then(|v| v.as_str())
        || src_meta.get("version").and_then(|v| v.as_i64()) != source_desc.get("version").and_then(|v| v.as_i64())
        || src_meta.get("profile").and_then(|v| v.as_str()) != source_ir.get("profile").and_then(|v| v.as_str())
        || src_meta.get("source_ir_sha256").and_then(|v| v.as_str())
            != source_desc.get("source_ir_sha256").and_then(|v| v.as_str())
        || src_meta.get("artifact_sha256").and_then(|v| v.as_str()) != Some(source_art_sha.as_str())
        || !src_meta
            .get("records")
            .unwrap()
            .exact_eq(source_ir.get("records").unwrap())
    {
        return Err("variant.source disagrees with source.json".to_string());
    }

    let variants_ir = payload.get("variants_ir").unwrap();
    expect_keys(variants_ir, &["variants"], "variants_ir")?;
    let variants = variants_ir
        .get("variants")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "variants_ir.variants must be an array".to_string())?;
    if variants.is_empty() || variants.len() > MAX_VARIANTS {
        return Err(format!("variants count must be 1..{MAX_VARIANTS}"));
    }

    let mut rec_by_id: HashMap<&str, &JsonValue> = HashMap::new();
    for rec in source_ir.get("records").unwrap().as_array().unwrap() {
        rec_by_id.insert(rec.get("record_id").unwrap().as_str().unwrap(), rec);
    }

    let mut seen_ids = HashSet::new();
    let mut allele_ids = Vec::with_capacity(variants.len());

    for (idx, entry) in variants.iter().enumerate() {
        let lbl = format!("variants[{idx}]");
        expect_keys(
            entry,
            &[
                "ordinal",
                "id",
                "record_id",
                "refget_id",
                "sequence_sha256",
                "edit_kind",
                "ref",
                "alt",
                "location",
                "state",
                "vrs_location_id",
                "vrs_allele_id",
                "vrs_allele",
                "mutated_bases",
                "mutated_sequence_sha256",
                "mutated_refget_id",
            ],
            &lbl,
        )?;
        if entry.get("ordinal").and_then(|v| v.as_i64()) != Some(idx as i64) {
            return Err(format!("{lbl}.ordinal mismatch"));
        }
        let vid = entry
            .get("id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("{lbl}.id must be string"))?;
        if vid.is_empty() || vid != vid.trim() || vid.len() > MAX_IDENTIFIER_BYTES {
            return Err(format!("{lbl}.id invalid"));
        }
        if !seen_ids.insert(vid.to_string()) {
            return Err(format!("duplicate variant id {vid:?}"));
        }

        let rid = entry
            .get("record_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("{lbl}.record_id must be string"))?;
        let rec = *rec_by_id
            .get(rid)
            .ok_or_else(|| format!("{lbl}.record_id {rid:?} absent from source"))?;
        let ref_seq = &sequences[rid];
        let rec_refget = rec.get("refget_id").unwrap().as_str().unwrap();
        let rec_sha = rec.get("sequence_sha256").unwrap().as_str().unwrap();
        if entry.get("refget_id").and_then(|v| v.as_str()) != Some(rec_refget)
            || entry.get("sequence_sha256").and_then(|v| v.as_str()) != Some(rec_sha)
        {
            return Err(format!("{lbl} refget_id/sequence_sha256 mismatch"));
        }

        let loc = entry.get("location").unwrap();
        expect_keys(loc, &["type", "sequenceReference", "start", "end"], &format!("{lbl}.location"))?;
        if loc.get("type").and_then(|v| v.as_str()) != Some("SequenceLocation") {
            return Err(format!("{lbl}.location.type invalid"));
        }
        let seq_ref = loc.get("sequenceReference").unwrap();
        expect_keys(seq_ref, &["type", "refgetAccession"], &format!("{lbl}.sequenceReference"))?;
        if seq_ref.get("type").and_then(|v| v.as_str()) != Some("SequenceReference")
            || seq_ref.get("refgetAccession").and_then(|v| v.as_str()) != Some(rec_refget)
        {
            return Err(format!("{lbl}.sequenceReference mismatch"));
        }
        let start = loc
            .get("start")
            .and_then(|v| v.as_i64())
            .ok_or_else(|| format!("{lbl}.start must be int"))?;
        let end = loc
            .get("end")
            .and_then(|v| v.as_i64())
            .ok_or_else(|| format!("{lbl}.end must be int"))?;
        if start < 0 || end < start || end as usize > ref_seq.len() {
            return Err(format!("{lbl} coordinates out of bounds"));
        }
        let (su, eu) = (start as usize, end as usize);

        let r_str = entry
            .get("ref")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("{lbl}.ref must be string"))?;
        let a_str = entry
            .get("alt")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("{lbl}.alt must be string"))?;
        if r_str.len() > MAX_VARIANT_BASES || a_str.len() > MAX_VARIANT_BASES {
            return Err(format!("{lbl} ref/alt exceeds base limit"));
        }
        if !r_str.bytes().all(is_iupac_base) || !a_str.bytes().all(is_iupac_base) {
            return Err(format!("{lbl} contains non-IUPAC bases"));
        }
        if su == eu && (!r_str.is_empty() || a_str.is_empty()) {
            return Err(format!("{lbl} insertion requires empty ref and nonempty alt"));
        }
        if su < eu && r_str.len() != eu - su {
            return Err(format!("{lbl} ref length disagrees with end - start"));
        }
        if &ref_seq[su..eu] != r_str {
            return Err(format!("{lbl} reference slice mismatch"));
        }
        if r_str == a_str {
            return Err(format!("{lbl} no-op edit forbidden"));
        }
        let expected_kind = match (r_str.len(), a_str.len()) {
            (0, x) if x > 0 => "insertion",
            (x, 0) if x > 0 => "deletion",
            (1, 1) => "snv",
            (x, y) if x == y && x > 1 => "mnv",
            _ => "delins",
        };
        if entry.get("edit_kind").and_then(|v| v.as_str()) != Some(expected_kind) {
            return Err(format!("{lbl}.edit_kind mismatch"));
        }

        let state = entry.get("state").unwrap();
        expect_keys(state, &["type", "sequence"], &format!("{lbl}.state"))?;
        if state.get("type").and_then(|v| v.as_str()) != Some("LiteralSequenceExpression")
            || state.get("sequence").and_then(|v| v.as_str()) != Some(a_str)
        {
            return Err(format!("{lbl}.state mismatch"));
        }

        let loc_dig = vrs_location_digest(rec_refget, start, end)?;
        let allele_dig = vrs_allele_digest(&loc_dig, a_str)?;
        let exp_loc_id = format!("ga4gh:SL.{loc_dig}");
        let exp_allele_id = format!("ga4gh:VA.{allele_dig}");
        let exp_vrs_allele = JsonValue::Object(vec![
            ("id".to_string(), JsonValue::String(exp_allele_id.clone())),
            ("type".to_string(), JsonValue::String("Allele".to_string())),
            ("digest".to_string(), JsonValue::String(allele_dig)),
            (
                "location".to_string(),
                JsonValue::Object(vec![
                    ("id".to_string(), JsonValue::String(exp_loc_id.clone())),
                    (
                        "type".to_string(),
                        JsonValue::String("SequenceLocation".to_string()),
                    ),
                    ("digest".to_string(), JsonValue::String(loc_dig)),
                    (
                        "sequenceReference".to_string(),
                        JsonValue::Object(vec![
                            (
                                "type".to_string(),
                                JsonValue::String("SequenceReference".to_string()),
                            ),
                            (
                                "refgetAccession".to_string(),
                                JsonValue::String(rec_refget.to_string()),
                            ),
                        ]),
                    ),
                    ("start".to_string(), JsonValue::Int(start)),
                    ("end".to_string(), JsonValue::Int(end)),
                ]),
            ),
            (
                "state".to_string(),
                JsonValue::Object(vec![
                    (
                        "type".to_string(),
                        JsonValue::String("LiteralSequenceExpression".to_string()),
                    ),
                    ("sequence".to_string(), JsonValue::String(a_str.to_string())),
                ]),
            ),
        ]);
        if entry.get("vrs_location_id").and_then(|v| v.as_str()) != Some(exp_loc_id.as_str())
            || entry.get("vrs_allele_id").and_then(|v| v.as_str()) != Some(exp_allele_id.as_str())
            || !entry.get("vrs_allele").unwrap().exact_eq(&exp_vrs_allele)
        {
            return Err(format!("{lbl} VRS 2.x identifier or object mismatch"));
        }

        let mut mutated = String::with_capacity(ref_seq.len() - (eu - su) + a_str.len());
        mutated.push_str(&ref_seq[..su]);
        mutated.push_str(a_str);
        mutated.push_str(&ref_seq[eu..]);
        let mut reverted = String::with_capacity(ref_seq.len());
        reverted.push_str(&mutated[..su]);
        reverted.push_str(r_str);
        reverted.push_str(&mutated[su + a_str.len()..]);
        if &reverted != ref_seq {
            return Err(format!("{lbl} reversion failed exact identity"));
        }
        let mut_sha = sha256_hex(mutated.as_bytes());
        let mut_refget = format!("SQ.{}", sha512t24u(mutated.as_bytes()));
        if entry.get("mutated_bases").and_then(|v| v.as_i64()) != Some(mutated.len() as i64)
            || entry.get("mutated_sequence_sha256").and_then(|v| v.as_str()) != Some(mut_sha.as_str())
            || entry.get("mutated_refget_id").and_then(|v| v.as_str()) != Some(mut_refget.as_str())
        {
            return Err(format!("{lbl} mutated sequence digest mismatch"));
        }
        allele_ids.push(JsonValue::String(exp_allele_id));
    }

    if payload.get("variants_ir_sha256").and_then(|v| v.as_str())
        != Some(digest_value(variants_ir)?.as_str())
    {
        return Err("variants_ir_sha256 mismatch".to_string());
    }
    let var_art_sha = verify_artifact_seal(&payload, "variant collection")?;

    let core = JsonValue::Object(vec![
        (
            "format".to_string(),
            JsonValue::String("brainc.variant-validation-report".to_string()),
        ),
        ("version".to_string(), JsonValue::Int(1)),
        (
            "validator".to_string(),
            JsonValue::Object(vec![
                (
                    "name".to_string(),
                    JsonValue::String("epigenesis-validate-rs".to_string()),
                ),
                ("version".to_string(), JsonValue::String("1.0.0".to_string())),
            ]),
        ),
        ("valid".to_string(), JsonValue::Bool(true)),
        ("code".to_string(), JsonValue::Null),
        ("detail".to_string(), JsonValue::Null),
        (
            "source_artifact_sha256".to_string(),
            JsonValue::String(source_art_sha),
        ),
        (
            "variant_artifact_sha256".to_string(),
            JsonValue::String(var_art_sha),
        ),
        (
            "checks".to_string(),
            JsonValue::Array(vec![
                JsonValue::Object(vec![
                    (
                        "name".to_string(),
                        JsonValue::String("source-bundle-integrity-and-sequences".to_string()),
                    ),
                    ("passed".to_string(), JsonValue::Bool(true)),
                ]),
                JsonValue::Object(vec![
                    (
                        "name".to_string(),
                        JsonValue::String("source-descriptor-binding".to_string()),
                    ),
                    ("passed".to_string(), JsonValue::Bool(true)),
                ]),
                JsonValue::Object(vec![
                    (
                        "name".to_string(),
                        JsonValue::String("vrs-v2-identifiers-and-reversion".to_string()),
                    ),
                    ("passed".to_string(), JsonValue::Bool(true)),
                ]),
                JsonValue::Object(vec![
                    (
                        "name".to_string(),
                        JsonValue::String("rfc8785-canonical-seals".to_string()),
                    ),
                    ("passed".to_string(), JsonValue::Bool(true)),
                ]),
            ]),
        ),
        (
            "summary".to_string(),
            JsonValue::Object(vec![
                ("variant_count".to_string(), JsonValue::Int(variants.len() as i64)),
                ("vrs_version".to_string(), JsonValue::String("2.x".to_string())),
                ("vrs_allele_ids".to_string(), JsonValue::Array(allele_ids)),
            ]),
        ),
    ]);
    let report_sha = digest_value(&core)?;
    let mut entries = core.as_object().unwrap().to_vec();
    entries.push(("report_sha256".to_string(), JsonValue::String(report_sha)));
    Ok(JsonValue::Object(entries))
}

pub fn validate_development_bundle_paths(
    source_dir: &Path,
    dev_dir: &Path,
    source_inputs: &BTreeMap<String, String>,
) -> Result<JsonValue, String> {
    let (source_desc, _children, _seqs) = verify_source_bundle_inner(source_dir, source_inputs)?;
    let source_art_sha = source_desc
        .get("artifact_sha256")
        .and_then(|v| v.as_str())
        .unwrap()
        .to_string();
    let source_ir_sha = source_desc
        .get("source_ir_sha256")
        .and_then(|v| v.as_str())
        .unwrap()
        .to_string();

    let children = [
        ("compilation_record", "compilation_record.json", "brainc.development-compilation", 1),
        ("development_module", "development_module.json", "brainc.development-module", 1),
        ("lowering_policy", "lowering_policy.json", "brainc.lowering-policy", 2),
        ("prediction_request", "prediction_request.json", "brainc.prediction-request", 2),
        ("prediction_response", "prediction_response.json", "brainc.prediction-response", 2),
        ("provider_manifest", "provider_manifest.json", "brainc.provider-manifest", 2),
        ("target_contract", "target_contract.json", "brainc.target-contract", 1),
    ];
    let mut limits = BTreeMap::new();
    limits.insert("bundle.json", MAX_JSON_BYTES);
    let mut required = vec!["bundle.json"];
    for &(_, fname, _, _) in &children {
        limits.insert(fname, MAX_JSON_BYTES);
        required.push(fname);
    }

    let raw_files = read_dir_files(dev_dir, &limits, &required, "development bundle")?;
    let bundle = parse_json_object(&raw_files["bundle.json"], "bundle.json")?;
    expect_keys(
        &bundle,
        &[
            "format",
            "version",
            "profile",
            "source",
            "artifacts",
            "blobs",
            "artifact_sha256",
        ],
        "bundle.json",
    )?;
    if bundle.get("format").and_then(|v| v.as_str()) != Some("brainc.development-bundle")
        || bundle.get("version").and_then(|v| v.as_i64()) != Some(1)
        || bundle.get("profile").and_then(|v| v.as_str()) != Some("development-module/v1")
    {
        return Err("bundle.json format/version/profile invalid".to_string());
    }
    let bundle_sha = verify_artifact_seal(&bundle, "bundle.json")?;

    let b_src = bundle
        .get("source")
        .ok_or_else(|| "bundle.json missing source".to_string())?;
    expect_keys(b_src, &["descriptor", "native"], "bundle.json.source")?;
    let b_src_desc = b_src.get("descriptor").unwrap();
    if b_src_desc.get("artifact_sha256").and_then(|v| v.as_str()) != Some(source_art_sha.as_str())
        || b_src_desc.get("ir_sha256").and_then(|v| v.as_str()) != Some(source_ir_sha.as_str())
    {
        return Err("bundle.json source.descriptor disagrees with source.json".to_string());
    }
    let source_ir_artifacts = source_desc
        .get("source_ir")
        .and_then(|s| s.get("artifacts"))
        .ok_or_else(|| "source.json missing source_ir.artifacts".to_string())?;
    if !b_src.get("native").unwrap().exact_eq(source_ir_artifacts) {
        return Err("bundle.json source.native disagrees with source.json".to_string());
    }

    let mut parsed_children = BTreeMap::new();
    let mut child_shas = BTreeMap::new();
    for &(role, fname, exp_fmt, exp_ver) in &children {
        let obj = parse_json_object(&raw_files[fname], fname)?;
        let sha = verify_artifact_seal(&obj, fname)?;
        if obj.get("format").and_then(|v| v.as_str()) != Some(exp_fmt)
            || obj.get("version").and_then(|v| v.as_i64()) != Some(exp_ver)
        {
            return Err(format!("{fname} format/version mismatch"));
        }
        let cref = bundle
            .get("artifacts")
            .and_then(|a| a.get(role))
            .ok_or_else(|| format!("bundle.json artifacts missing {role}"))?;
        if cref.get("format").and_then(|v| v.as_str()) != Some(exp_fmt)
            || cref.get("version").and_then(|v| v.as_i64()) != Some(exp_ver)
            || cref.get("artifact_sha256").and_then(|v| v.as_str()) != Some(sha.as_str())
        {
            return Err(format!("bundle.json artifacts.{role} disagrees with {fname}"));
        }
        parsed_children.insert(role, obj);
        child_shas.insert(role, sha);
    }

    // Verify development_module internal module_sha256, tensors, budgets, and blobs closure
    let dev_mod = &parsed_children["development_module"];
    let module_inner = dev_mod
        .get("module")
        .ok_or_else(|| "development_module.json missing module".to_string())?;
    let module_sha = digest_value(module_inner)?;
    if dev_mod.get("module_sha256").and_then(|v| v.as_str()) != Some(module_sha.as_str()) {
        return Err("development_module.module_sha256 mismatch".to_string());
    }
    let tensors = module_inner
        .get("tensors")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "development_module.module.tensors must be array".to_string())?;
    let mut expected_blob_shas = BTreeSet::new();
    let mut total_tensor_bytes: i64 = 0;
    for (idx, t) in tensors.iter().enumerate() {
        let storage = t
            .get("storage")
            .ok_or_else(|| format!("tensor[{idx}] missing storage"))?;
        let kind = storage
            .get("kind")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("tensor[{idx}] storage missing kind"))?;
        let exp_len = storage
            .get("byte_length")
            .and_then(|v| v.as_i64())
            .ok_or_else(|| format!("tensor[{idx}] missing byte_length"))?;
        let exp_sha = storage
            .get("sha256")
            .and_then(|v| v.as_str())
            .ok_or_else(|| format!("tensor[{idx}] missing sha256"))?;
        total_tensor_bytes += exp_len;
        if kind == "inline-base64" {
            let b64 = storage
                .get("data")
                .and_then(|v| v.as_str())
                .ok_or_else(|| format!("tensor[{idx}] inline-base64 missing data"))?;
            let decoded = base64_decode_strict(b64)?;
            if decoded.len() as i64 != exp_len || sha256_hex(&decoded) != exp_sha {
                return Err(format!("tensor[{idx}] inline-base64 digest/length mismatch"));
            }
        } else if kind == "sha256-blob" {
            if !is_sha256_hex(exp_sha) {
                return Err(format!("tensor[{idx}] invalid sha256-blob digest"));
            }
            expected_blob_shas.insert(exp_sha.to_string());
            return Err("external sha256-blob verification requires external blob store".to_string());
        } else {
            return Err(format!("tensor[{idx}] unsupported storage kind {kind:?}"));
        }
    }

    // Verify bundle.json blobs closure
    let blobs_arr = bundle
        .get("blobs")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "bundle.json blobs must be an array".to_string())?;
    if !blobs_arr.is_empty() || !expected_blob_shas.is_empty() {
        return Err("bundle.json blobs mismatch or unverified external blobs".to_string());
    }

    // Verify prediction_response outputs and inline-base64 digests
    let man = &parsed_children["provider_manifest"];
    let req = &parsed_children["prediction_request"];
    let resp = &parsed_children["prediction_response"];
    let pol = &parsed_children["lowering_policy"];
    let tgt = &parsed_children["target_contract"];

    let resp_outputs = resp
        .get("outputs")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "prediction_response missing outputs array".to_string())?;
    for (idx, out) in resp_outputs.iter().enumerate() {
        let storage = out
            .get("storage")
            .ok_or_else(|| format!("prediction_response.outputs[{idx}] missing storage"))?;
        if storage.get("kind").and_then(|v| v.as_str()) == Some("inline-base64") {
            let b64 = storage
                .get("data")
                .and_then(|v| v.as_str())
                .ok_or_else(|| format!("response.outputs[{idx}] missing data"))?;
            let decoded = base64_decode_strict(b64)?;
            let exp_len = storage
                .get("byte_length")
                .and_then(|v| v.as_i64())
                .ok_or_else(|| format!("response.outputs[{idx}] missing byte_length"))?;
            let exp_sha = storage
                .get("sha256")
                .and_then(|v| v.as_str())
                .ok_or_else(|| format!("response.outputs[{idx}] missing sha256"))?;
            if decoded.len() as i64 != exp_len || sha256_hex(&decoded) != exp_sha {
                return Err(format!("response.outputs[{idx}] inline-base64 digest/length mismatch"));
            }
        }
    }

    // Verify manifest/request/response cross-references and source acceptance
    let src_profile = source_desc
        .get("source_ir")
        .and_then(|s| s.get("profile"))
        .and_then(|v| v.as_str())
        .unwrap_or("");
    let accepts = man
        .get("accepts")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "provider_manifest missing accepts array".to_string())?;
    let expected_accept_prefix = format!("brainc.source-descriptor/v1;profile={src_profile}");
    if !accepts
        .iter()
        .any(|a| a.as_str().map_or(false, |s| s == expected_accept_prefix || s.starts_with(&format!("{expected_accept_prefix};"))))
    {
        return Err("provider_manifest.accepts does not accept source profile".to_string());
    }
    if req.get("provider_manifest_sha256").and_then(|v| v.as_str())
        != Some(child_shas["provider_manifest"].as_str())
        || req.get("source").and_then(|s| s.get("artifact_sha256")).and_then(|v| v.as_str())
            != Some(source_art_sha.as_str())
        || resp.get("request_artifact_sha256").and_then(|v| v.as_str())
            != Some(child_shas["prediction_request"].as_str())
        || !resp
            .get("provider")
            .unwrap_or(&JsonValue::Null)
            .exact_eq(man.get("provider").unwrap_or(&JsonValue::Null))
        || !resp
            .get("model_identity")
            .unwrap_or(&JsonValue::Null)
            .exact_eq(man.get("model_identity").unwrap_or(&JsonValue::Null))
    {
        return Err("prediction request/response/manifest binding mismatch".to_string());
    }

    // Verify lowering_policy -> target_contract binding and module budgets/compilation_record
    let tgt_contract = tgt
        .get("contract")
        .ok_or_else(|| "target_contract missing contract".to_string())?;
    let tgt_contract_sha = digest_value(tgt_contract)?;
    if tgt.get("contract_sha256").and_then(|v| v.as_str()) != Some(tgt_contract_sha.as_str()) {
        return Err("target_contract.contract_sha256 mismatch".to_string());
    }
    if pol
        .get("target")
        .and_then(|t| t.get("contract_sha256"))
        .and_then(|v| v.as_str())
        != Some(tgt_contract_sha.as_str())
    {
        return Err("lowering_policy target.contract_sha256 disagrees with target_contract".to_string());
    }

    let ops = module_inner
        .get("entrypoint")
        .and_then(|ep| ep.get("operations"))
        .and_then(|v| v.as_array())
        .ok_or_else(|| "development_module.module.entrypoint.operations must be array".to_string())?;
    let budgets = module_inner
        .get("budgets")
        .ok_or_else(|| "development_module.module missing budgets".to_string())?;
    if budgets.get("operations").and_then(|v| v.as_i64()) != Some(ops.len() as i64)
        || budgets.get("tensor_bytes").and_then(|v| v.as_i64()) != Some(total_tensor_bytes)
    {
        return Err("development_module.module.budgets mismatch".to_string());
    }

    // Verify compilation_record cross-references and budgets
    let comp = &parsed_children["compilation_record"];
    if comp.get("source").and_then(|s| s.get("artifact_sha256")).and_then(|v| v.as_str())
        != Some(source_art_sha.as_str())
    {
        return Err("compilation_record.source disagrees with source.json".to_string());
    }
    let comp_artifacts = comp
        .get("artifacts")
        .ok_or_else(|| "compilation_record missing artifacts".to_string())?;
    for &(role, _, exp_fmt, exp_ver) in &children {
        if role == "compilation_record" {
            continue;
        }
        let cref = comp_artifacts
            .get(role)
            .ok_or_else(|| format!("compilation_record.artifacts missing {role}"))?;
        if cref.get("format").and_then(|v| v.as_str()) != Some(exp_fmt)
            || cref.get("version").and_then(|v| v.as_i64()) != Some(exp_ver)
            || cref.get("artifact_sha256").and_then(|v| v.as_str()) != Some(child_shas[role].as_str())
        {
            return Err(format!("compilation_record.artifacts.{role} mismatch"));
        }
    }
    if comp
        .get("result")
        .and_then(|r| r.get("module_sha256"))
        .and_then(|v| v.as_str())
        != Some(module_sha.as_str())
        || !comp
            .get("result")
            .and_then(|r| r.get("budgets"))
            .unwrap_or(&JsonValue::Null)
            .exact_eq(budgets)
    {
        return Err("compilation_record.result disagrees with development_module".to_string());
    }

    let _ = jcs_string;
    Ok(JsonValue::Object(vec![
        ("valid".to_string(), JsonValue::Bool(true)),
        (
            "validation_scope".to_string(),
            JsonValue::String("source-closure+structural-and-seal-replay".to_string()),
        ),
        (
            "validator".to_string(),
            JsonValue::String("epigenesis-validate-rs/1.0.0".to_string()),
        ),
        (
            "source_artifact_sha256".to_string(),
            JsonValue::String(source_art_sha),
        ),
        (
            "development_bundle_artifact_sha256".to_string(),
            JsonValue::String(bundle_sha),
        ),
        (
            "development_module_artifact_sha256".to_string(),
            JsonValue::String(child_shas["development_module"].clone()),
        ),
        ("module_sha256".to_string(), JsonValue::String(module_sha)),
        (
            "compilation_artifact_sha256".to_string(),
            JsonValue::String(child_shas["compilation_record"].clone()),
        ),
    ]))
}
