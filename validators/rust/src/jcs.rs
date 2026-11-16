//! Strict I-JSON parser and RFC 8785 (JCS) canonicalizer.

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::crypto::{sha256_hex, sha512t24u};

pub const SAFE_INTEGER: i64 = (1_i64 << 53) - 1;
pub const MAX_JSON_DEPTH: usize = 64;
pub const MAX_JSON_MEMBERS: usize = 1_000_000;
pub const MAX_STRING_BYTES: usize = 1024 * 1024;

#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    String(String),
    Array(Vec<JsonValue>),
    Object(Vec<(String, JsonValue)>),
}

impl JsonValue {
    pub fn as_object(&self) -> Option<&[(String, JsonValue)]> {
        match self {
            JsonValue::Object(entries) => Some(entries),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[JsonValue]> {
        match self {
            JsonValue::Array(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            JsonValue::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            JsonValue::Int(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            JsonValue::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        self.as_object()
            .and_then(|entries| entries.iter().find(|(k, _)| k == key).map(|(_, v)| v))
    }

    pub fn exact_eq(&self, other: &Self) -> bool {
        match (self, other) {
            (JsonValue::Null, JsonValue::Null) => true,
            (JsonValue::Bool(a), JsonValue::Bool(b)) => a == b,
            (JsonValue::Int(a), JsonValue::Int(b)) => a == b,
            (JsonValue::Float(a), JsonValue::Float(b)) => a.to_bits() == b.to_bits() || (*a == 0.0 && *b == 0.0),
            (JsonValue::String(a), JsonValue::String(b)) => a == b,
            (JsonValue::Array(a), JsonValue::Array(b)) => {
                a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| x.exact_eq(y))
            }
            (JsonValue::Object(a), JsonValue::Object(b)) => {
                if a.len() != b.len() {
                    return false;
                }
                let mut ma: BTreeMap<&str, &JsonValue> = BTreeMap::new();
                for (k, v) in a {
                    ma.insert(k.as_str(), v);
                }
                let mut mb: BTreeMap<&str, &JsonValue> = BTreeMap::new();
                for (k, v) in b {
                    mb.insert(k.as_str(), v);
                }
                ma.len() == mb.len()
                    && ma.iter().all(|(k, va)| mb.get(k).is_some_and(|vb| va.exact_eq(vb)))
            }
            _ => false,
        }
    }
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
    members: usize,
    max_members: usize,
}

impl<'a> Parser<'a> {
    fn new(bytes: &'a [u8], max_members: usize) -> Self {
        Self {
            bytes,
            pos: 0,
            members: 0,
            max_members,
        }
    }

    fn skip_ws(&mut self) {
        while self.pos < self.bytes.len() {
            match self.bytes[self.pos] {
                b' ' | b'\t' | b'\r' | b'\n' => self.pos += 1,
                _ => break,
            }
        }
    }

    fn bump(&mut self) -> Result<u8, String> {
        if self.pos >= self.bytes.len() {
            return Err("unexpected end of JSON input".to_string());
        }
        let b = self.bytes[self.pos];
        self.pos += 1;
        Ok(b)
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn count_member(&mut self) -> Result<(), String> {
        self.members += 1;
        if self.members > self.max_members {
            return Err("JSON exceeds maximum member count".to_string());
        }
        Ok(())
    }

    fn parse_value(&mut self, depth: usize) -> Result<JsonValue, String> {
        if depth > MAX_JSON_DEPTH {
            return Err("JSON exceeds maximum nesting depth".to_string());
        }
        self.count_member()?;
        self.skip_ws();
        let b = self.peek().ok_or_else(|| "unexpected end of JSON".to_string())?;
        match b {
            b'{' => self.parse_object(depth),
            b'[' => self.parse_array(depth),
            b'"' => self.parse_string().map(JsonValue::String),
            b't' => self.parse_literal(b"true", JsonValue::Bool(true)),
            b'f' => self.parse_literal(b"false", JsonValue::Bool(false)),
            b'n' => self.parse_literal(b"null", JsonValue::Null),
            b'-' | b'0'..=b'9' => self.parse_number(),
            _ => Err(format!("unexpected byte 0x{b:02x} in JSON")),
        }
    }

    fn parse_literal(&mut self, lit: &[u8], val: JsonValue) -> Result<JsonValue, String> {
        if self.bytes[self.pos..].starts_with(lit) {
            self.pos += lit.len();
            Ok(val)
        } else {
            Err("invalid JSON literal".to_string())
        }
    }

    fn parse_hex4(&mut self) -> Result<u16, String> {
        let mut val: u16 = 0;
        for _ in 0..4 {
            let b = self.bump()?;
            let nibble = match b {
                b'0'..=b'9' => (b - b'0') as u16,
                b'a'..=b'f' => (b - b'a' + 10) as u16,
                b'A'..=b'F' => (b - b'A' + 10) as u16,
                _ => return Err("invalid hex digit in \\u escape".to_string()),
            };
            val = (val << 4) | nibble;
        }
        Ok(val)
    }

    fn parse_string(&mut self) -> Result<String, String> {
        if self.bump()? != b'"' {
            return Err("expected '\"'".to_string());
        }
        let mut out = String::new();
        while self.pos < self.bytes.len() {
            let b = self.bytes[self.pos];
            if b == b'"' {
                self.pos += 1;
                if out.len() > MAX_STRING_BYTES {
                    return Err("JSON string exceeds maximum byte limit".to_string());
                }
                return Ok(out);
            }
            if b < 0x20 {
                return Err("unescaped control character in JSON string".to_string());
            }
            if b == b'\\' {
                self.pos += 1;
                let esc = self.bump()?;
                match esc {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'/' => out.push('/'),
                    b'b' => out.push('\x08'),
                    b'f' => out.push('\x0c'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        let code = self.parse_hex4()?;
                        if (0xD800..=0xDBFF).contains(&code) {
                            if self.pos + 2 > self.bytes.len()
                                || self.bytes[self.pos] != b'\\'
                                || self.bytes[self.pos + 1] != b'u'
                            {
                                return Err("lone high Unicode surrogate in JSON string".to_string());
                            }
                            self.pos += 2;
                            let low = self.parse_hex4()?;
                            if !(0xDC00..=0xDFFF).contains(&low) {
                                return Err("invalid low Unicode surrogate in JSON string".to_string());
                            }
                            let cp = 0x10000 + (((code as u32 - 0xD800) << 10) | (low as u32 - 0xDC00));
                            let ch = char::from_u32(cp)
                                .ok_or_else(|| "invalid Unicode scalar value".to_string())?;
                            out.push(ch);
                        } else if (0xDC00..=0xDFFF).contains(&code) {
                            return Err("lone low Unicode surrogate in JSON string".to_string());
                        } else {
                            let ch = char::from_u32(code as u32)
                                .ok_or_else(|| "invalid Unicode scalar value".to_string())?;
                            out.push(ch);
                        }
                    }
                    _ => return Err("invalid escape sequence in JSON string".to_string()),
                }
            } else {
                let start = self.pos;
                while self.pos < self.bytes.len() {
                    let b = self.bytes[self.pos];
                    if b == b'"' || b == b'\\' || b < 0x20 {
                        break;
                    }
                    self.pos += 1;
                }
                let slice = std::str::from_utf8(&self.bytes[start..self.pos])
                    .map_err(|_| "invalid UTF-8 in JSON string".to_string())?;
                out.push_str(slice);
            }
            if out.len() > MAX_STRING_BYTES {
                return Err("JSON string exceeds maximum byte limit".to_string());
            }
        }
        Err("unterminated JSON string".to_string())
    }

    fn parse_number(&mut self) -> Result<JsonValue, String> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek() {
            Some(b'0') => {
                self.pos += 1;
                if matches!(self.peek(), Some(b'0'..=b'9')) {
                    return Err("leading zero in JSON number".to_string());
                }
            }
            Some(b'1'..=b'9') => {
                self.pos += 1;
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.pos += 1;
                }
            }
            _ => return Err("invalid JSON number".to_string()),
        }
        let mut is_float = false;
        if self.peek() == Some(b'.') {
            is_float = true;
            self.pos += 1;
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err("expected fraction digits in JSON number".to_string());
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            is_float = true;
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err("expected exponent digits in JSON number".to_string());
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        let slice = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| "invalid UTF-8 in number".to_string())?;
        if is_float {
            let val: f64 = slice.parse().map_err(|_| "invalid float literal".to_string())?;
            if !val.is_finite() {
                return Err("non-finite float in JSON".to_string());
            }
            Ok(JsonValue::Float(val))
        } else {
            let val: i64 = slice
                .parse()
                .map_err(|_| format!("integer outside safe range: {slice}"))?;
            if !(-SAFE_INTEGER..=SAFE_INTEGER).contains(&val) {
                return Err(format!("integer outside safe range: {val}"));
            }
            Ok(JsonValue::Int(val))
        }
    }

    fn parse_array(&mut self, depth: usize) -> Result<JsonValue, String> {
        self.bump()?; // '['
        self.skip_ws();
        let mut items = Vec::new();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(JsonValue::Array(items));
        }
        loop {
            let item = self.parse_value(depth + 1)?;
            items.push(item);
            self.skip_ws();
            match self.bump()? {
                b',' => continue,
                b']' => return Ok(JsonValue::Array(items)),
                _ => return Err("expected ',' or ']' in JSON array".to_string()),
            }
        }
    }

    fn parse_object(&mut self, depth: usize) -> Result<JsonValue, String> {
        self.bump()?; // '{'
        self.skip_ws();
        let mut entries = Vec::new();
        let mut seen = std::collections::HashSet::new();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(JsonValue::Object(entries));
        }
        loop {
            self.skip_ws();
            if self.peek() != Some(b'"') {
                return Err("JSON object key must be a string".to_string());
            }
            let key = self.parse_string()?;
            if !seen.insert(key.clone()) {
                return Err(format!("duplicate JSON key: {key:?}"));
            }
            self.skip_ws();
            if self.bump()? != b':' {
                return Err("expected ':' after JSON object key".to_string());
            }
            let val = self.parse_value(depth + 1)?;
            entries.push((key, val));
            self.skip_ws();
            match self.bump()? {
                b',' => continue,
                b'}' => return Ok(JsonValue::Object(entries)),
                _ => return Err("expected ',' or '}' in JSON object".to_string()),
            }
        }
    }
}

pub fn parse_json_object(raw: &[u8], label: &str) -> Result<JsonValue, String> {
    parse_json_object_with_limit(raw, label, MAX_JSON_MEMBERS)
}

pub fn parse_json_object_with_limit(
    raw: &[u8],
    label: &str,
    max_members: usize,
) -> Result<JsonValue, String> {
    let text = std::str::from_utf8(raw).map_err(|e| format!("{label} is not valid UTF-8: {e}"))?;
    let mut parser = Parser::new(text.as_bytes(), max_members);
    let value = parser
        .parse_value(0)
        .map_err(|e| format!("malformed JSON in {label}: {e}"))?;
    parser.skip_ws();
    if parser.pos != parser.bytes.len() {
        return Err(format!("trailing content after JSON in {label}"));
    }
    if !matches!(value, JsonValue::Object(_)) {
        return Err(format!("{label} must be a top-level JSON object"));
    }
    Ok(value)
}

pub fn jcs_string(s: &str, out: &mut String) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\x08' => out.push_str("\\b"),
            '\x0c' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

pub fn jcs_number(value: f64, out: &mut String) -> Result<(), String> {
    if !value.is_finite() {
        return Err("non-finite number in JSON".to_string());
    }
    if value == 0.0 {
        out.push('0');
        return Ok(());
    }
    let abs = value.abs();
    if value < 0.0 {
        out.push('-');
    }
    // Rust {:e} produces exact shortest roundtrip decimal significand and exponent.
    let raw = format!("{abs:e}");
    let (mantissa, exp_str) = raw.split_once('e').expect("scientific format contains e");
    let exponent: i32 = exp_str.parse().expect("valid exponent");
    let dot_idx = mantissa.find('.');
    let mut digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let int_len = dot_idx.unwrap_or(mantissa.len()) as i32;
    let mut point = int_len + exponent;
    while digits.len() > 1 && digits.starts_with('0') {
        digits.remove(0);
        point -= 1;
    }
    if (1e-6..1e21).contains(&abs) {
        let mut rendered = String::new();
        if point <= 0 {
            rendered.push_str("0.");
            for _ in 0..(-point) {
                rendered.push('0');
            }
            rendered.push_str(&digits);
        } else if (point as usize) >= digits.len() {
            rendered.push_str(&digits);
            for _ in 0..(point as usize - digits.len()) {
                rendered.push('0');
            }
        } else {
            let p = point as usize;
            rendered.push_str(&digits[..p]);
            rendered.push('.');
            rendered.push_str(&digits[p..]);
        }
        if rendered.contains('.') {
            while rendered.ends_with('0') {
                rendered.pop();
            }
            if rendered.ends_with('.') {
                rendered.pop();
            }
        }
        out.push_str(&rendered);
        return Ok(());
    }
    while digits.ends_with('0') && digits.len() > 1 {
        digits.pop();
    }
    let norm_exp = point - 1;
    out.push(digits.as_bytes()[0] as char);
    if digits.len() > 1 {
        out.push('.');
        out.push_str(&digits[1..]);
    }
    out.push('e');
    if norm_exp >= 0 {
        out.push('+');
    }
    let _ = write!(out, "{norm_exp}");
    Ok(())
}

fn canonical_write(value: &JsonValue, depth: usize, out: &mut String) -> Result<(), String> {
    if depth > MAX_JSON_DEPTH {
        return Err("JSON exceeds maximum depth during canonicalization".to_string());
    }
    match value {
        JsonValue::Null => out.push_str("null"),
        JsonValue::Bool(true) => out.push_str("true"),
        JsonValue::Bool(false) => out.push_str("false"),
        JsonValue::Int(n) => {
            if !(-SAFE_INTEGER..=SAFE_INTEGER).contains(n) {
                return Err(format!("integer outside safe range: {n}"));
            }
            let _ = write!(out, "{n}");
        }
        JsonValue::Float(v) => jcs_number(*v, out)?,
        JsonValue::String(s) => jcs_string(s, out),
        JsonValue::Array(items) => {
            out.push('[');
            for (idx, item) in items.iter().enumerate() {
                if idx > 0 {
                    out.push(',');
                }
                canonical_write(item, depth + 1, out)?;
            }
            out.push(']');
        }
        JsonValue::Object(entries) => {
            let mut sorted: Vec<(&(String, JsonValue), Vec<u16>)> = entries
                .iter()
                .map(|pair| (pair, pair.0.encode_utf16().collect::<Vec<u16>>()))
                .collect();
            sorted.sort_by(|a, b| a.1.cmp(&b.1));
            out.push('{');
            for (idx, ((k, v), _)) in sorted.iter().enumerate() {
                if idx > 0 {
                    out.push(',');
                }
                jcs_string(k, out);
                out.push(':');
                canonical_write(v, depth + 1, out)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

pub fn canonical_bytes(value: &JsonValue) -> Result<Vec<u8>, String> {
    let mut out = String::new();
    canonical_write(value, 0, &mut out)?;
    Ok(out.into_bytes())
}

pub fn digest_value(value: &JsonValue) -> Result<String, String> {
    Ok(sha256_hex(&canonical_bytes(value)?))
}

pub fn digest_object_without(entries: &[(String, JsonValue)], omitted_key: &str) -> Result<String, String> {
    let filtered: Vec<(String, JsonValue)> = entries
        .iter()
        .filter(|(k, _)| k != omitted_key)
        .cloned()
        .collect();
    digest_value(&JsonValue::Object(filtered))
}

pub fn sha512t24u_value(value: &JsonValue) -> Result<String, String> {
    Ok(sha512t24u(&canonical_bytes(value)?))
}
