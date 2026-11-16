pub mod crypto;
pub mod gzip;
pub mod jcs;
pub mod validator;

#[cfg(test)]
mod tests {
    use super::crypto::{sha256_hex, sha512t24u};
    use super::jcs::{canonical_bytes, parse_json_object};
    use super::validator::{vrs_allele_digest, vrs_location_digest};

    #[test]
    fn sha256_known_vector() {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn ga4gh_vrs_v2_apoe_vector() {
        let loc = vrs_location_digest("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 44908821, 44908822)
            .expect("valid location digest");
        assert_eq!(loc, "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz");
        let allele = vrs_allele_digest(&loc, "T").expect("valid allele digest");
        assert_eq!(allele, "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
    }

    #[test]
    fn rfc8785_sorting_and_duplicate_rejection() {
        let parsed = parse_json_object(
            br#"{"b":1.0,"a":-0.0,"u\ud83d\ude00":1e-6,"\u00e9":1e21}"#,
            "test",
        )
        .expect("valid JSON");
        let canon = String::from_utf8(canonical_bytes(&parsed).expect("canonical")).unwrap();
        assert_eq!(canon, r#"{"a":0,"b":1,"u😀":0.000001,"é":1e+21}"#);

        assert!(parse_json_object(br#"{"k":1,"k":2}"#, "dup").is_err());
        assert!(parse_json_object(br#"{"k":NaN}"#, "nan").is_err());
        assert!(parse_json_object(br#"{"k":"\ud800"}"#, "surrogate").is_err());
        assert!(parse_json_object(br#"{"k":9007199254740992}"#, "unsafe_int").is_err());
    }

    #[test]
    fn sha512t24u_sequence_id() {
        assert_eq!(
            format!("SQ.{}", sha512t24u(b"ACGT")),
            "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
        );
    }
}
