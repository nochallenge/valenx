//! RON envelope for round-trip persistence of a [`crate::Schematic`].

use serde::{Deserialize, Serialize};

use crate::error::SymbolError;
use crate::schematic::Schematic;

/// File format version. Bumped on schema changes.
pub const VERSION: u32 = 1;

/// On-wire envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SchematicFile {
    /// Format version.
    pub version: u32,
    /// Schematic payload.
    pub schematic: Schematic,
}

impl SchematicFile {
    /// Wrap a schematic in the current envelope.
    pub fn new(schematic: Schematic) -> Self {
        Self {
            version: VERSION,
            schematic,
        }
    }
}

/// Serialise a schematic to a pretty RON string.
pub fn to_ron_string(s: &Schematic) -> Result<String, SymbolError> {
    let file = SchematicFile::new(s.clone());
    ron::ser::to_string_pretty(&file, ron::ser::PrettyConfig::default())
        .map_err(|e| SymbolError::Ron(e.to_string()))
}

/// Parse a schematic from a RON string. Fails on version mismatch
/// rather than silently coercing.
pub fn from_ron_str(s: &str) -> Result<Schematic, SymbolError> {
    let file: SchematicFile = ron::de::from_str(s).map_err(|e| SymbolError::Ron(e.to_string()))?;
    if file.version != VERSION {
        return Err(SymbolError::Ron(format!(
            "version mismatch: file = {}, expected = {}",
            file.version, VERSION
        )));
    }
    Ok(file.schematic)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schematic::{PlacedSymbol, Wire};
    use crate::symbol::SymbolKind;

    #[test]
    fn round_trips_small_schematic() {
        let mut s = Schematic::new();
        s.push_symbol(PlacedSymbol::new(SymbolKind::Resistor, [0.0, 0.0]));
        s.push_wire(Wire::new(vec![[0.0, 0.0], [10.0, 0.0]], "N1").unwrap());
        let text = to_ron_string(&s).expect("ser");
        let back = from_ron_str(&text).expect("de");
        assert_eq!(back.symbols.len(), 1);
        assert_eq!(back.wires.len(), 1);
    }

    #[test]
    fn rejects_version_mismatch() {
        let bad = "(version: 99, schematic: (symbols: [], wires: []))";
        assert!(matches!(from_ron_str(bad), Err(SymbolError::Ron(_))));
    }
}
