//! Catalogs — indexed views over `Field` and `ScalarRecord`.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::field::Field;
use crate::scalar::ScalarRecord;
use crate::time::TimeKey;

/// Indexed collection of `Field` values. Fields are keyed by name
/// and (implicitly) by their `TimeKey`. Lookups by name over a time
/// series are the common case.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct FieldCatalog {
    entries: Vec<Field>,
}

impl FieldCatalog {
    /// New, empty catalog.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Append a field (without de-duplicating by `(name, time)`).
    pub fn insert(&mut self, field: Field) {
        self.entries.push(field);
    }

    /// Total entry count (one per stored field).
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// `true` when the catalog has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate distinct field names in insertion order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        // De-duplicate names in insertion order.
        let mut seen = std::collections::HashSet::new();
        self.entries
            .iter()
            .filter_map(move |f| seen.insert(f.name.clone()).then_some(f.name.as_str()))
    }

    /// All entries with the given name (typically a time series).
    pub fn by_name<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Field> + 'a {
        self.entries.iter().filter(move |f| f.name == name)
    }

    /// Look up one field by `(name, time)`. Returns the first match.
    pub fn at_time(&self, name: &str, time: TimeKey) -> Option<&Field> {
        self.entries
            .iter()
            .find(|f| f.name == name && f.time == time)
    }

    /// All `TimeKey`s for which a named field has data, in insertion
    /// order.
    pub fn time_series(&self, name: &str) -> Vec<TimeKey> {
        self.entries
            .iter()
            .filter(|f| f.name == name)
            .map(|f| f.time)
            .collect()
    }
}

/// Indexed collection of scalar records, keyed by name.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ScalarCatalog {
    entries: BTreeMap<String, Vec<ScalarRecord>>,
}

impl ScalarCatalog {
    /// New, empty catalog.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a scalar record under its name (appending to the
    /// per-name `Vec` when one already exists).
    pub fn insert(&mut self, record: ScalarRecord) {
        self.entries
            .entry(record.name.clone())
            .or_default()
            .push(record);
    }

    /// First record (in insertion order) registered under `name`, or
    /// `None` if the name is unknown.
    pub fn get(&self, name: &str) -> Option<&ScalarRecord> {
        self.entries.get(name).and_then(|v| v.first())
    }

    /// All records registered under `name`, in insertion order; empty
    /// slice when the name is unknown.
    pub fn all(&self, name: &str) -> &[ScalarRecord] {
        self.entries.get(name).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Iterate the unique names present in the catalog (lexical order).
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.keys().map(|k| k.as_str())
    }

    /// Total record count summed across all names.
    pub fn len(&self) -> usize {
        self.entries.values().map(|v| v.len()).sum()
    }

    /// `true` when no records have been inserted.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Field, FieldKind, Location, RegionRef};
    use crate::scalar::{ScalarRecord, ScalarSource};
    use crate::time::TimeKey;
    use crate::units::{DIMENSIONLESS, PASCAL};

    fn make_scalar_field(name: &str, time: TimeKey) -> Field {
        Field {
            name: name.into(),
            kind: FieldKind::Scalar,
            location: Location::OnCell,
            region: RegionRef("fluid".into()),
            units: PASCAL,
            time,
            data: vec![],
            range: None,
        }
    }

    #[test]
    fn field_catalog_time_series() {
        let mut cat = FieldCatalog::new();
        cat.insert(make_scalar_field("p", TimeKey::Iteration(1)));
        cat.insert(make_scalar_field("p", TimeKey::Iteration(2)));
        cat.insert(make_scalar_field("u", TimeKey::Steady));
        assert_eq!(cat.time_series("p").len(), 2);
        assert_eq!(cat.time_series("u").len(), 1);
        assert!(cat.at_time("p", TimeKey::Iteration(2)).is_some());
    }

    #[test]
    fn scalar_catalog_insert_and_get() {
        let mut cat = ScalarCatalog::new();
        cat.insert(ScalarRecord {
            name: "cd".into(),
            value: 0.031,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Computed,
            description: None,
        });
        assert_eq!(cat.get("cd").unwrap().value, 0.031);
        assert!(cat.get("cl").is_none());
    }
}
