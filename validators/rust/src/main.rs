use std::collections::BTreeMap;
use std::env;
use std::path::Path;
use std::process;

use epigenesis_validate_rs::jcs::{canonical_bytes, digest_value, parse_json_object, sha512t24u_value};
use epigenesis_validate_rs::validator::{
    read_regular_file, validate_development_bundle_paths, validate_variant_collection_paths,
    MAX_JSON_BYTES,
};

fn parse_source_inputs(args: &[String]) -> Result<(Vec<String>, BTreeMap<String, String>), String> {
    let mut positional = Vec::new();
    let mut source_inputs = BTreeMap::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--source-input" {
            i += 1;
            let spec = args
                .get(i)
                .ok_or_else(|| "--source-input requires ROLE=PATH".to_string())?;
            let (role, path) = spec
                .split_once('=')
                .ok_or_else(|| "--source-input must be ROLE=PATH".to_string())?;
            if source_inputs.insert(role.to_string(), path.to_string()).is_some() {
                return Err(format!("duplicate --source-input role {role:?}"));
            }
        } else {
            positional.push(args[i].clone());
        }
        i += 1;
    }
    Ok((positional, source_inputs))
}

fn run() -> Result<i32, String> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        return Err("usage: epigenesis-validate-rs <validate-variants|validate-development|canonical-digest> ...".to_string());
    }
    match args[1].as_str() {
        "canonical-digest" => {
            if args.len() != 3 {
                return Err("usage: epigenesis-validate-rs canonical-digest <file.json>".to_string());
            }
            let raw = read_regular_file(Path::new(&args[2]), "json input", MAX_JSON_BYTES)?;
            let val = parse_json_object(&raw, "json input")?;
            let sha256 = digest_value(&val)?;
            let t24u = sha512t24u_value(&val)?;
            println!(r#"{{"sha256":"{sha256}","sha512t24u":"{t24u}"}}"#);
            Ok(0)
        }
        "validate-variants" => {
            let (pos, inputs) = parse_source_inputs(&args[2..])?;
            if pos.len() != 2 {
                return Err("usage: epigenesis-validate-rs validate-variants <source_dir> <variant_collection.json> [--source-input ROLE=PATH]".to_string());
            }
            match validate_variant_collection_paths(Path::new(&pos[0]), Path::new(&pos[1]), &inputs) {
                Ok(report) => {
                    let bytes = canonical_bytes(&report)?;
                    println!("{}", String::from_utf8(bytes).unwrap());
                    Ok(0)
                }
                Err(report_json) => {
                    println!("{report_json}");
                    Ok(3)
                }
            }
        }
        "validate-development" => {
            let (pos, inputs) = parse_source_inputs(&args[2..])?;
            if pos.len() != 2 {
                return Err("usage: epigenesis-validate-rs validate-development <source_dir> <development_dir> [--source-input ROLE=PATH]".to_string());
            }
            let report =
                validate_development_bundle_paths(Path::new(&pos[0]), Path::new(&pos[1]), &inputs)?;
            let bytes = canonical_bytes(&report)?;
            println!("{}", String::from_utf8(bytes).unwrap());
            Ok(0)
        }
        cmd => Err(format!("unknown subcommand: {cmd}")),
    }
}

fn main() {
    match run() {
        Ok(code) => process::exit(code),
        Err(msg) => {
            eprintln!("epigenesis-validate-rs: {msg}");
            process::exit(2);
        }
    }
}
