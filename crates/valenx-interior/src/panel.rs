//! UI panel envelope — rooms + furniture palette + add-by-click placeholder.

use nalgebra::Vector3;

use crate::error::InteriorError;
use crate::furniture::{Furniture, Placement};
use crate::room::Room;

/// Workbench panel state.
#[derive(Default)]
pub struct InteriorPanelState {
    /// All rooms in the scene.
    pub rooms: Vec<Room>,
    /// All placed furniture instances.
    pub placements: Vec<Placement>,
    /// Currently-selected palette item.
    pub palette_selection: Option<Furniture>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl InteriorPanelState {
    /// New empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a room — fails on duplicate id.
    pub fn add_room(&mut self, room: Room) -> Result<(), InteriorError> {
        if self.rooms.iter().any(|r| r.id == room.id) {
            return Err(InteriorError::BadParameter {
                name: "room.id",
                reason: format!("duplicate id `{}`", room.id),
            });
        }
        self.last_status = Some(format!("added room `{}`", room.id));
        self.last_error = None;
        self.rooms.push(room);
        Ok(())
    }

    /// Select a palette item — what the next click-to-place spawns.
    pub fn select(&mut self, f: Furniture) {
        self.palette_selection = Some(f);
        self.last_status = Some(format!("selected `{}`", f.name()));
        self.last_error = None;
    }

    /// Add-by-click placeholder — spawns the selected item at the
    /// requested world position inside the named room.
    pub fn click_to_place(
        &mut self,
        position: Vector3<f64>,
        room_id: &str,
    ) -> Result<(), InteriorError> {
        let kind = self.palette_selection.ok_or(InteriorError::BadParameter {
            name: "palette",
            reason: "no palette selection".into(),
        })?;
        if !self.rooms.iter().any(|r| r.id == room_id) {
            return Err(InteriorError::UnknownRoom(room_id.into()));
        }
        let placement = kind.place(position, 0.0, room_id);
        self.last_status = Some(format!("placed `{}` in `{}`", kind.name(), room_id));
        self.last_error = None;
        self.placements.push(placement);
        Ok(())
    }

    /// Status setter.
    pub fn set_status(&mut self, s: impl Into<String>) {
        self.last_status = Some(s.into());
        self.last_error = None;
    }

    /// Error setter.
    pub fn set_error(&mut self, s: impl Into<String>) {
        self.last_error = Some(s.into());
        self.last_status = None;
    }
}
