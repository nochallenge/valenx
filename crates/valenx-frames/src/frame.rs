//! A frame = a collection of members + their joints.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::FramesError;
use crate::member::Member;

/// A joint reference — which member endpoints share this joint.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Joint {
    /// World position of the joint (mm).
    pub position: Vector3<f64>,
    /// (member_index, endpoint) — `endpoint` is 0 for the first
    /// path vertex of the member, len(path)-1 for the last.
    pub connected: Vec<(usize, usize)>,
}

/// A structural frame.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Frame {
    /// All members.
    pub members: Vec<Member>,
    /// Joints inferred from coincident member endpoints (or
    /// supplied explicitly).
    pub joints: Vec<Joint>,
}

impl Frame {
    /// Empty frame.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a member, returning its index.
    pub fn push_member(&mut self, m: Member) -> usize {
        self.members.push(m);
        self.members.len() - 1
    }

    /// Append an explicit joint, returning its index.
    pub fn push_joint(&mut self, j: Joint) -> Result<usize, FramesError> {
        for &(mid, _) in &j.connected {
            if mid >= self.members.len() {
                return Err(FramesError::BadIndex {
                    got: mid,
                    n: self.members.len(),
                });
            }
        }
        self.joints.push(j);
        Ok(self.joints.len() - 1)
    }

    /// Total cut-list length (mm).
    pub fn total_length_mm(&self) -> f64 {
        self.members.iter().map(|m| m.length_mm()).sum()
    }

    /// Auto-detect joints from coincident endpoints within `tol`
    /// (mm). Overwrites the joints vector.
    pub fn auto_joints(&mut self, tol: f64) {
        self.joints.clear();
        let tol2 = tol * tol;
        let mut endpoints: Vec<(usize, usize, Vector3<f64>)> = Vec::new();
        for (mi, m) in self.members.iter().enumerate() {
            let last = m.path.len() - 1;
            endpoints.push((mi, 0, m.path[0]));
            endpoints.push((mi, last, *m.path.last().unwrap()));
        }
        let mut consumed = vec![false; endpoints.len()];
        for i in 0..endpoints.len() {
            if consumed[i] {
                continue;
            }
            let (mi, ep, p) = endpoints[i];
            let mut conn = vec![(mi, ep)];
            consumed[i] = true;
            for j in (i + 1)..endpoints.len() {
                if consumed[j] {
                    continue;
                }
                let (mj, ej, q) = endpoints[j];
                if (p - q).norm_squared() <= tol2 {
                    conn.push((mj, ej));
                    consumed[j] = true;
                }
            }
            if conn.len() > 1 {
                self.joints.push(Joint {
                    position: p,
                    connected: conn,
                });
            }
        }
    }
}

/// RON envelope version.
pub const VERSION: u32 = 1;

/// Persistence envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FrameFile {
    /// Format version.
    pub version: u32,
    /// Frame payload.
    pub frame: Frame,
}

/// Serialise a frame to RON.
pub fn to_ron_string(f: &Frame) -> Result<String, FramesError> {
    let file = FrameFile {
        version: VERSION,
        frame: f.clone(),
    };
    ron::ser::to_string_pretty(&file, ron::ser::PrettyConfig::default())
        .map_err(|e| FramesError::Ron(e.to_string()))
}

/// Parse a frame from RON.
pub fn from_ron_str(s: &str) -> Result<Frame, FramesError> {
    let file: FrameFile = ron::de::from_str(s).map_err(|e| FramesError::Ron(e.to_string()))?;
    if file.version != VERSION {
        return Err(FramesError::Ron(format!(
            "version mismatch: file = {}, expected = {}",
            file.version, VERSION
        )));
    }
    Ok(file.frame)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile::Profile;

    #[test]
    fn auto_joints_detects_l_intersection() {
        let mut f = Frame::new();
        f.push_member(Member::straight(
            Vector3::zeros(),
            Vector3::new(1000.0, 0.0, 0.0),
            Profile::default_ipe200(),
        ));
        f.push_member(Member::straight(
            Vector3::new(1000.0, 0.0, 0.0),
            Vector3::new(1000.0, 1000.0, 0.0),
            Profile::default_ipe200(),
        ));
        f.auto_joints(1.0);
        assert_eq!(f.joints.len(), 1);
        assert_eq!(f.joints[0].connected.len(), 2);
    }

    #[test]
    fn total_length_matches_cut_list() {
        let mut f = Frame::new();
        f.push_member(Member::straight(
            Vector3::zeros(),
            Vector3::new(2000.0, 0.0, 0.0),
            Profile::default_ipe200(),
        ));
        f.push_member(Member::straight(
            Vector3::new(2000.0, 0.0, 0.0),
            Vector3::new(2000.0, 1500.0, 0.0),
            Profile::default_ipe200(),
        ));
        assert!((f.total_length_mm() - 3500.0).abs() < 1e-6);
    }

    #[test]
    fn ron_round_trips() {
        let mut f = Frame::new();
        f.push_member(Member::straight(
            Vector3::zeros(),
            Vector3::new(1000.0, 0.0, 0.0),
            Profile::default_ipe200(),
        ));
        let text = to_ron_string(&f).unwrap();
        let back = from_ron_str(&text).unwrap();
        assert_eq!(back.members.len(), 1);
    }
}
