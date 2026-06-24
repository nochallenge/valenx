//! IEEE 1278.1 Distributed Interactive Simulation (DIS) — **Entity State
//! PDU** codec.
//!
//! Honest scope: this implements *only* the Entity State PDU (PDU type 1),
//! which is by far the most-used PDU and the natural one for a co-simulation
//! federate to publish/consume. It is **not** the full PDU family. The byte
//! layout is bit-exact to the standard: all multi-byte fields are
//! **big-endian** (network byte order), the fixed record is **144 bytes**
//! (with zero articulation parameters), and each articulation parameter
//! record adds 16 bytes.
//!
//! Field layout (offsets are from the start of the PDU):
//!
//! ```text
//! PDU header (12 bytes)
//!   0   u8   protocol version          (6 = IEEE 1278.1-1995, 7 = -2012)
//!   1   u8   exercise id
//!   2   u8   pdu type                  (1 = Entity State)
//!   3   u8   protocol family           (1 = Entity Information/Interaction)
//!   4   u32  timestamp
//!   8   u16  pdu length (bytes)
//!  10   u8   pdu status
//!  11   u8   padding
//! Entity State body
//!  12   EntityId (6)  site u16, application u16, entity u16
//!  18   u8   force id
//!  19   u8   number of articulation parameters (N)
//!  20   EntityType (8)  kind u8, domain u8, country u16,
//!                       category u8, subcategory u8, specific u8, extra u8
//!  28   EntityType (8)  alternative entity type, same layout
//!  36   Vec3f32 (12)  linear velocity x,y,z
//!  48   Vec3f64 (24)  location x,y,z   (double precision)
//!  72   Vec3f32 (12)  orientation psi,theta,phi (radians)
//!  84   u32  appearance
//!  88   DeadReckoning (40)  algorithm u8, 15 bytes other params,
//!                           linear accel 3xf32, angular velocity 3xf32
//! 128   EntityMarking (12)  charset u8, 11 marking bytes
//! 140   u32  capabilities
//! 144   articulation parameters: N x 16 bytes
//! ```

use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use std::io::{Cursor, Write};

use crate::error::{FmiError, Result};

/// PDU type number for the Entity State PDU.
pub const PDU_TYPE_ENTITY_STATE: u8 = 1;
/// Protocol family for Entity Information / Interaction.
pub const PROTOCOL_FAMILY_ENTITY_INFO: u8 = 1;
/// Byte length of the Entity State PDU with zero articulation parameters.
pub const ENTITY_STATE_FIXED_LEN: usize = 144;
/// Byte length of one articulation (variable) parameter record.
pub const ARTICULATION_PARAM_LEN: usize = 16;

/// DIS Entity Identifier: `(site, application, entity)`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EntityId {
    /// Simulation site number.
    pub site: u16,
    /// Application number within the site.
    pub application: u16,
    /// Entity number within the application.
    pub entity: u16,
}

/// DIS Entity Type (a.k.a. entity type record), the standard 7-field
/// taxonomy that classifies an entity.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EntityType {
    /// Kind (platform, munition, lifeform, …).
    pub kind: u8,
    /// Domain (land, air, surface, subsurface, space).
    pub domain: u8,
    /// Country code.
    pub country: u16,
    /// Category within kind/domain.
    pub category: u8,
    /// Subcategory.
    pub subcategory: u8,
    /// Specific.
    pub specific: u8,
    /// Extra.
    pub extra: u8,
}

/// A three-component single-precision vector (velocity / orientation /
/// acceleration / angular velocity).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3f32 {
    /// X component.
    pub x: f32,
    /// Y component.
    pub y: f32,
    /// Z component.
    pub z: f32,
}

/// A three-component double-precision vector (world location).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3f64 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

/// Dead-reckoning parameters record (40 bytes).
///
/// The all-zero default corresponds to dead-reckoning algorithm 0 ("Other")
/// with no acceleration or angular velocity — a valid, inert record.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DeadReckoning {
    /// Dead-reckoning algorithm number.
    pub algorithm: u8,
    /// 15 bytes of "other parameters" (quaternion / unused per algorithm).
    pub other_parameters: [u8; 15],
    /// Entity linear acceleration.
    pub linear_acceleration: Vec3f32,
    /// Entity angular velocity.
    pub angular_velocity: Vec3f32,
}

/// Entity marking record (12 bytes): a character-set byte plus 11 bytes of
/// marking (typically the entity's name, space-padded).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EntityMarking {
    /// Character set identifier (1 = ASCII).
    pub character_set: u8,
    /// 11 marking characters.
    pub marking: [u8; 11],
}

impl Default for EntityMarking {
    fn default() -> Self {
        Self {
            character_set: 1,
            marking: [0u8; 11],
        }
    }
}

impl EntityMarking {
    /// Build an ASCII marking record from `s`, truncating or
    /// space-padding to 11 bytes.
    pub fn ascii(s: &str) -> Self {
        let mut marking = [0u8; 11];
        for (slot, b) in marking.iter_mut().zip(s.bytes()) {
            *slot = b;
        }
        Self {
            character_set: 1,
            marking,
        }
    }
}

/// One articulation (variable) parameter record (16 bytes). Opaque here —
/// the bytes are carried verbatim so a round-trip is exact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ArticulationParameter {
    /// The raw 16-byte record.
    pub bytes: [u8; ARTICULATION_PARAM_LEN],
}

/// A fully-decoded DIS Entity State PDU.
#[derive(Clone, Debug, PartialEq)]
pub struct EntityStatePdu {
    /// Protocol version (6 or 7).
    pub protocol_version: u8,
    /// Exercise identifier.
    pub exercise_id: u8,
    /// Timestamp field (units per the standard's timestamp scheme).
    pub timestamp: u32,
    /// PDU status byte.
    pub pdu_status: u8,
    /// Entity identifier.
    pub entity_id: EntityId,
    /// Force identifier (friendly / opposing / neutral / other).
    pub force_id: u8,
    /// Entity type taxonomy.
    pub entity_type: EntityType,
    /// Alternative (perceived) entity type.
    pub alternative_entity_type: EntityType,
    /// Entity linear velocity (m/s).
    pub linear_velocity: Vec3f32,
    /// Entity location in world coordinates (m).
    pub location: Vec3f64,
    /// Entity orientation Euler angles (rad).
    pub orientation: Vec3f32,
    /// Appearance bitfield.
    pub appearance: u32,
    /// Dead-reckoning parameters.
    pub dead_reckoning: DeadReckoning,
    /// Entity marking.
    pub marking: EntityMarking,
    /// Capabilities bitfield.
    pub capabilities: u32,
    /// Articulation parameters (each 16 bytes).
    pub articulation_parameters: Vec<ArticulationParameter>,
}

impl Default for EntityStatePdu {
    fn default() -> Self {
        Self {
            protocol_version: 6,
            exercise_id: 0,
            timestamp: 0,
            pdu_status: 0,
            entity_id: EntityId::default(),
            force_id: 0,
            entity_type: EntityType::default(),
            alternative_entity_type: EntityType::default(),
            linear_velocity: Vec3f32::default(),
            location: Vec3f64::default(),
            orientation: Vec3f32::default(),
            appearance: 0,
            dead_reckoning: DeadReckoning::default(),
            marking: EntityMarking::default(),
            capabilities: 0,
            articulation_parameters: Vec::new(),
        }
    }
}

impl EntityStatePdu {
    /// Total serialized length in bytes (fixed 144 + 16 per articulation
    /// parameter).
    pub fn byte_len(&self) -> usize {
        ENTITY_STATE_FIXED_LEN + self.articulation_parameters.len() * ARTICULATION_PARAM_LEN
    }

    /// Serialize this PDU to its exact IEEE 1278.1 big-endian byte layout.
    ///
    /// The PDU header's `pdu type`, `protocol family`, and `pdu length`
    /// fields are written authoritatively from the structure (length =
    /// [`EntityStatePdu::byte_len`], type = [`PDU_TYPE_ENTITY_STATE`]),
    /// and the articulation-parameter count byte is written from the
    /// actual number of records, so a serialized PDU is always internally
    /// consistent.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.byte_len());

        // --- PDU header (12 bytes) ---
        // Infallible: writing into a Vec never errors. `expect` documents
        // that invariant rather than silently ignoring a Result.
        buf.write_u8(self.protocol_version).expect("vec write");
        buf.write_u8(self.exercise_id).expect("vec write");
        buf.write_u8(PDU_TYPE_ENTITY_STATE).expect("vec write");
        buf.write_u8(PROTOCOL_FAMILY_ENTITY_INFO)
            .expect("vec write");
        buf.write_u32::<BigEndian>(self.timestamp)
            .expect("vec write");
        buf.write_u16::<BigEndian>(self.byte_len() as u16)
            .expect("vec write");
        buf.write_u8(self.pdu_status).expect("vec write");
        buf.write_u8(0).expect("vec write"); // header padding

        // --- Entity ID (6 bytes) ---
        buf.write_u16::<BigEndian>(self.entity_id.site)
            .expect("vec write");
        buf.write_u16::<BigEndian>(self.entity_id.application)
            .expect("vec write");
        buf.write_u16::<BigEndian>(self.entity_id.entity)
            .expect("vec write");

        // --- Force ID + articulation parameter count (2 bytes) ---
        buf.write_u8(self.force_id).expect("vec write");
        buf.write_u8(self.articulation_parameters.len() as u8)
            .expect("vec write");

        // --- Entity type + alternative entity type (8 + 8 bytes) ---
        write_entity_type(&mut buf, &self.entity_type);
        write_entity_type(&mut buf, &self.alternative_entity_type);

        // --- Linear velocity (12 bytes) ---
        write_vec3f32(&mut buf, &self.linear_velocity);

        // --- Location (24 bytes, double precision) ---
        write_vec3f64(&mut buf, &self.location);

        // --- Orientation (12 bytes) ---
        write_vec3f32(&mut buf, &self.orientation);

        // --- Appearance (4 bytes) ---
        buf.write_u32::<BigEndian>(self.appearance)
            .expect("vec write");

        // --- Dead reckoning (40 bytes) ---
        buf.write_u8(self.dead_reckoning.algorithm)
            .expect("vec write");
        buf.write_all(&self.dead_reckoning.other_parameters)
            .expect("vec write");
        write_vec3f32(&mut buf, &self.dead_reckoning.linear_acceleration);
        write_vec3f32(&mut buf, &self.dead_reckoning.angular_velocity);

        // --- Entity marking (12 bytes) ---
        buf.write_u8(self.marking.character_set).expect("vec write");
        buf.write_all(&self.marking.marking).expect("vec write");

        // --- Capabilities (4 bytes) ---
        buf.write_u32::<BigEndian>(self.capabilities)
            .expect("vec write");

        // --- Articulation parameters (16 bytes each) ---
        for ap in &self.articulation_parameters {
            buf.write_all(&ap.bytes).expect("vec write");
        }

        debug_assert_eq!(buf.len(), self.byte_len());
        buf
    }

    /// Decode an Entity State PDU from its exact byte layout.
    ///
    /// Fail-loud: a buffer shorter than the fixed record, a header that
    /// declares a non-Entity-State PDU type, or an articulation-parameter
    /// count that overruns the buffer all return [`FmiError`] rather than
    /// reading garbage or panicking.
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < ENTITY_STATE_FIXED_LEN {
            return Err(FmiError::PduTooShort {
                needed: ENTITY_STATE_FIXED_LEN,
                got: data.len(),
                what: "Entity State PDU fixed record",
            });
        }
        let mut c = Cursor::new(data);

        // --- PDU header ---
        let protocol_version = read_u8(&mut c)?;
        let exercise_id = read_u8(&mut c)?;
        let pdu_type = read_u8(&mut c)?;
        if pdu_type != PDU_TYPE_ENTITY_STATE {
            return Err(FmiError::UnsupportedPduType {
                got: pdu_type,
                expected: PDU_TYPE_ENTITY_STATE,
            });
        }
        let _protocol_family = read_u8(&mut c)?;
        let timestamp = read_u32(&mut c)?;
        let _pdu_length = read_u16(&mut c)?;
        let pdu_status = read_u8(&mut c)?;
        let _padding = read_u8(&mut c)?;

        // --- Entity ID ---
        let entity_id = EntityId {
            site: read_u16(&mut c)?,
            application: read_u16(&mut c)?,
            entity: read_u16(&mut c)?,
        };

        // --- Force ID + articulation parameter count ---
        let force_id = read_u8(&mut c)?;
        let n_articulation = read_u8(&mut c)? as usize;

        // --- Entity type records ---
        let entity_type = read_entity_type(&mut c)?;
        let alternative_entity_type = read_entity_type(&mut c)?;

        // --- Linear velocity ---
        let linear_velocity = read_vec3f32(&mut c)?;

        // --- Location ---
        let location = read_vec3f64(&mut c)?;

        // --- Orientation ---
        let orientation = read_vec3f32(&mut c)?;

        // --- Appearance ---
        let appearance = read_u32(&mut c)?;

        // --- Dead reckoning ---
        let algorithm = read_u8(&mut c)?;
        let mut other_parameters = [0u8; 15];
        read_exact(&mut c, &mut other_parameters)?;
        let linear_acceleration = read_vec3f32(&mut c)?;
        let angular_velocity = read_vec3f32(&mut c)?;
        let dead_reckoning = DeadReckoning {
            algorithm,
            other_parameters,
            linear_acceleration,
            angular_velocity,
        };

        // --- Entity marking ---
        let character_set = read_u8(&mut c)?;
        let mut marking_bytes = [0u8; 11];
        read_exact(&mut c, &mut marking_bytes)?;
        let marking = EntityMarking {
            character_set,
            marking: marking_bytes,
        };

        // --- Capabilities ---
        let capabilities = read_u32(&mut c)?;

        // --- Articulation parameters ---
        let needed = ENTITY_STATE_FIXED_LEN + n_articulation * ARTICULATION_PARAM_LEN;
        if data.len() < needed {
            return Err(FmiError::PduTooShort {
                needed,
                got: data.len(),
                what: "articulation parameter records",
            });
        }
        let mut articulation_parameters = Vec::with_capacity(n_articulation);
        for _ in 0..n_articulation {
            let mut bytes = [0u8; ARTICULATION_PARAM_LEN];
            read_exact(&mut c, &mut bytes)?;
            articulation_parameters.push(ArticulationParameter { bytes });
        }

        Ok(EntityStatePdu {
            protocol_version,
            exercise_id,
            timestamp,
            pdu_status,
            entity_id,
            force_id,
            entity_type,
            alternative_entity_type,
            linear_velocity,
            location,
            orientation,
            appearance,
            dead_reckoning,
            marking,
            capabilities,
            articulation_parameters,
        })
    }
}

// --- write helpers (big-endian) ---

fn write_entity_type(buf: &mut Vec<u8>, t: &EntityType) {
    buf.write_u8(t.kind).expect("vec write");
    buf.write_u8(t.domain).expect("vec write");
    buf.write_u16::<BigEndian>(t.country).expect("vec write");
    buf.write_u8(t.category).expect("vec write");
    buf.write_u8(t.subcategory).expect("vec write");
    buf.write_u8(t.specific).expect("vec write");
    buf.write_u8(t.extra).expect("vec write");
}

fn write_vec3f32(buf: &mut Vec<u8>, v: &Vec3f32) {
    buf.write_f32::<BigEndian>(v.x).expect("vec write");
    buf.write_f32::<BigEndian>(v.y).expect("vec write");
    buf.write_f32::<BigEndian>(v.z).expect("vec write");
}

fn write_vec3f64(buf: &mut Vec<u8>, v: &Vec3f64) {
    buf.write_f64::<BigEndian>(v.x).expect("vec write");
    buf.write_f64::<BigEndian>(v.y).expect("vec write");
    buf.write_f64::<BigEndian>(v.z).expect("vec write");
}

// --- read helpers (big-endian, fail-loud on short reads) ---

fn short(what: &'static str) -> FmiError {
    FmiError::PduTooShort {
        needed: 0,
        got: 0,
        what,
    }
}

fn read_u8(c: &mut Cursor<&[u8]>) -> Result<u8> {
    c.read_u8().map_err(|_| short("u8 field"))
}
fn read_u16(c: &mut Cursor<&[u8]>) -> Result<u16> {
    c.read_u16::<BigEndian>().map_err(|_| short("u16 field"))
}
fn read_u32(c: &mut Cursor<&[u8]>) -> Result<u32> {
    c.read_u32::<BigEndian>().map_err(|_| short("u32 field"))
}
fn read_f32(c: &mut Cursor<&[u8]>) -> Result<f32> {
    c.read_f32::<BigEndian>().map_err(|_| short("f32 field"))
}
fn read_f64(c: &mut Cursor<&[u8]>) -> Result<f64> {
    c.read_f64::<BigEndian>().map_err(|_| short("f64 field"))
}
fn read_exact(c: &mut Cursor<&[u8]>, out: &mut [u8]) -> Result<()> {
    use std::io::Read;
    c.read_exact(out).map_err(|_| short("byte block"))
}

fn read_entity_type(c: &mut Cursor<&[u8]>) -> Result<EntityType> {
    Ok(EntityType {
        kind: read_u8(c)?,
        domain: read_u8(c)?,
        country: read_u16(c)?,
        category: read_u8(c)?,
        subcategory: read_u8(c)?,
        specific: read_u8(c)?,
        extra: read_u8(c)?,
    })
}

fn read_vec3f32(c: &mut Cursor<&[u8]>) -> Result<Vec3f32> {
    Ok(Vec3f32 {
        x: read_f32(c)?,
        y: read_f32(c)?,
        z: read_f32(c)?,
    })
}

fn read_vec3f64(c: &mut Cursor<&[u8]>) -> Result<Vec3f64> {
    Ok(Vec3f64 {
        x: read_f64(c)?,
        y: read_f64(c)?,
        z: read_f64(c)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_pdu() -> EntityStatePdu {
        EntityStatePdu {
            protocol_version: 6,
            exercise_id: 7,
            timestamp: 0x0102_0304,
            pdu_status: 0,
            entity_id: EntityId {
                site: 1,
                application: 2,
                entity: 3,
            },
            force_id: 1,
            entity_type: EntityType {
                kind: 1,
                domain: 1,
                country: 225,
                category: 1,
                subcategory: 1,
                specific: 0,
                extra: 0,
            },
            alternative_entity_type: EntityType::default(),
            linear_velocity: Vec3f32 {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            },
            location: Vec3f64 {
                x: 100.5,
                y: -200.25,
                z: 300.125,
            },
            orientation: Vec3f32 {
                x: 0.1,
                y: 0.2,
                z: 0.3,
            },
            appearance: 0xDEAD_BEEF,
            dead_reckoning: DeadReckoning {
                algorithm: 2,
                other_parameters: [0u8; 15],
                linear_acceleration: Vec3f32 {
                    x: 0.5,
                    y: 0.0,
                    z: -9.81,
                },
                angular_velocity: Vec3f32::default(),
            },
            marking: EntityMarking::ascii("TANK01"),
            capabilities: 0,
            articulation_parameters: Vec::new(),
        }
    }

    #[test]
    fn fixed_pdu_is_144_bytes() {
        let pdu = sample_pdu();
        assert_eq!(pdu.byte_len(), 144);
        assert_eq!(pdu.to_bytes().len(), 144);
    }

    #[test]
    fn round_trips_exactly() {
        let pdu = sample_pdu();
        let bytes = pdu.to_bytes();
        let decoded = EntityStatePdu::from_bytes(&bytes).expect("decode");
        assert_eq!(decoded, pdu);
    }

    #[test]
    fn round_trips_with_articulation_parameters() {
        let mut pdu = sample_pdu();
        pdu.articulation_parameters
            .push(ArticulationParameter { bytes: [0xAB; 16] });
        pdu.articulation_parameters
            .push(ArticulationParameter { bytes: [0xCD; 16] });
        assert_eq!(pdu.byte_len(), 144 + 32);
        let bytes = pdu.to_bytes();
        assert_eq!(bytes.len(), 176);
        let decoded = EntityStatePdu::from_bytes(&bytes).expect("decode");
        assert_eq!(decoded, pdu);
    }

    /// Pin the EXACT bytes of specific header / body offsets against the
    /// IEEE 1278.1 layout. This is the standard-conformance anchor: if the
    /// layout ever shifts, these byte assertions fail loudly.
    #[test]
    fn reference_pdu_byte_layout_is_exact() {
        let pdu = sample_pdu();
        let b = pdu.to_bytes();

        // Header.
        assert_eq!(b[0], 6, "protocol version @0");
        assert_eq!(b[1], 7, "exercise id @1");
        assert_eq!(b[2], PDU_TYPE_ENTITY_STATE, "pdu type @2 = 1");
        assert_eq!(b[3], PROTOCOL_FAMILY_ENTITY_INFO, "protocol family @3 = 1");
        // Timestamp 0x01020304 big-endian @4..8.
        assert_eq!(&b[4..8], &[0x01, 0x02, 0x03, 0x04], "timestamp @4 BE");
        // PDU length 144 = 0x0090 big-endian @8..10.
        assert_eq!(&b[8..10], &[0x00, 0x90], "pdu length @8 BE = 144");
        assert_eq!(b[10], 0, "pdu status @10");
        assert_eq!(b[11], 0, "header padding @11");

        // Entity ID: site 1, application 2, entity 3 — each u16 BE.
        assert_eq!(&b[12..14], &[0x00, 0x01], "site @12 BE = 1");
        assert_eq!(&b[14..16], &[0x00, 0x02], "application @14 BE = 2");
        assert_eq!(&b[16..18], &[0x00, 0x03], "entity @16 BE = 3");

        // Force ID @18, articulation count @19 (zero here).
        assert_eq!(b[18], 1, "force id @18");
        assert_eq!(b[19], 0, "articulation count @19 = 0");

        // Entity type @20: kind 1, domain 1, country 225 (0x00E1) BE.
        assert_eq!(b[20], 1, "entity kind @20");
        assert_eq!(b[21], 1, "entity domain @21");
        assert_eq!(&b[22..24], &[0x00, 0xE1], "country @22 BE = 225");

        // Linear velocity @36: x=1.0,y=2.0,z=3.0 as BE f32.
        assert_eq!(&b[36..40], &1.0f32.to_be_bytes(), "vel.x @36 BE f32");
        assert_eq!(&b[40..44], &2.0f32.to_be_bytes(), "vel.y @40 BE f32");
        assert_eq!(&b[44..48], &3.0f32.to_be_bytes(), "vel.z @44 BE f32");

        // Location @48: double precision BE f64.
        assert_eq!(&b[48..56], &100.5f64.to_be_bytes(), "loc.x @48 BE f64");
        assert_eq!(&b[56..64], &(-200.25f64).to_be_bytes(), "loc.y @56 BE f64");
        assert_eq!(&b[64..72], &300.125f64.to_be_bytes(), "loc.z @64 BE f64");

        // Orientation @72: BE f32.
        assert_eq!(&b[72..76], &0.1f32.to_be_bytes(), "ori.psi @72 BE f32");

        // Appearance @84: 0xDEADBEEF BE.
        assert_eq!(&b[84..88], &[0xDE, 0xAD, 0xBE, 0xEF], "appearance @84 BE");

        // Dead reckoning algorithm @88.
        assert_eq!(b[88], 2, "dead-reckoning algorithm @88");

        // Entity marking @128: charset 1 @128, "TANK01" @129..
        assert_eq!(b[128], 1, "marking charset @128 = ASCII");
        assert_eq!(&b[129..135], b"TANK01", "marking text @129");
        assert_eq!(&b[135..140], &[0u8; 5], "marking pad @135");

        // Capabilities @140..144.
        assert_eq!(&b[140..144], &[0x00, 0x00, 0x00, 0x00], "capabilities @140");
    }

    #[test]
    fn short_buffer_is_fail_loud() {
        let bytes = vec![0u8; 100]; // < 144
        assert!(matches!(
            EntityStatePdu::from_bytes(&bytes),
            Err(FmiError::PduTooShort { .. })
        ));
    }

    #[test]
    fn wrong_pdu_type_is_fail_loud() {
        let mut bytes = sample_pdu().to_bytes();
        bytes[2] = 2; // claim "Fire" PDU, not Entity State
        assert!(matches!(
            EntityStatePdu::from_bytes(&bytes),
            Err(FmiError::UnsupportedPduType {
                got: 2,
                expected: 1
            })
        ));
    }

    #[test]
    fn articulation_count_overrun_is_fail_loud() {
        let mut bytes = sample_pdu().to_bytes();
        // Claim 5 articulation parameters but supply none of the 80 bytes.
        bytes[19] = 5;
        assert!(matches!(
            EntityStatePdu::from_bytes(&bytes),
            Err(FmiError::PduTooShort { .. })
        ));
    }
}
