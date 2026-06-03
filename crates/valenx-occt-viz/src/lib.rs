//! # valenx-occt-viz
//!
//! Round 5 Block 4 (Phases 161-200) — OpenCASCADE (OCCT) **visualization
//! patterns** feature parity for Valenx. **This crate closes the entire
//! 200-phase FreeCAD-parity roadmap.**
//!
//! OCCT's visualization subsystem ships as the V3d / AIS / Prs3d / SelectMgr
//! cluster: roughly 100k LOC of C++ that drives the OpenGl rendering
//! pipeline plus a HUD-style object selection and presentation model
//! ([interactive context], [highlight states], [Drawer aspect bundles]).
//!
//! Direct-porting that layer is the wrong move for Valenx — the existing
//! [`valenx-app`] shell already runs an `egui` + `wgpu` viewport with
//! [`OrbitCamera`], depth-buffer-based hidden-surface removal, and
//! per-frame mesh upload via [`WgpuRenderer`]. This crate **adapts** the
//! OCCT V3d / AIS / Prs3d concepts onto that pipeline rather than
//! shipping a second renderer alongside it.
//!
//! - **`V3d_View_*`** camera + viewer ops (orbit, pan, zoom, fit-all,
//!   axonometric snaps, perspective↔orthographic, depth buffer,
//!   back-face culling, clipping planes, XOR drag overlay).
//! - **`AIS_InteractiveContext` + `SelectMgr`** selection ops (single,
//!   box, polygon, hover/selected highlight, face/edge/vertex picking,
//!   subshape-kind filters).
//! - **`Prs3d_Drawer`** display-aspect bundles (material presets,
//!   per-face/edge colour, transparency, line width / style, hidden-line
//!   display, isoparametric U/V lines, section-plane outlines).
//! - **`AIS_Animation_*` + `gp_Trsf`** animation + transformation
//!   manipulators (camera path tween, object pose interpolation,
//!   exploded view, XYZ axis drag widget, rotation-ring widget).
//! - **Misc viz** screenshot capture, video-frame export, navigation
//!   cube widget, scalar-field legend, world-origin axes marker, XY
//!   ground-plane grid.
//!
//! This crate provides a Rust-native API surface for 40 of those, one
//! module per OCCT-equivalent feature so each can advance from
//! "scaffold" → "v1 implementation" → "production parity"
//! independently.
//!
//! [interactive context]: https://dev.opencascade.org/doc/refman/html/class_a_i_s___interactive_context.html
//! [highlight states]: https://dev.opencascade.org/doc/refman/html/class_a_i_s___interactive_object.html
//! [Drawer aspect bundles]: https://dev.opencascade.org/doc/refman/html/class_prs3d___drawer.html
//! [`valenx-app`]: ../valenx_app/index.html
//! [`OrbitCamera`]: ../valenx_viz/camera/struct.OrbitCamera.html
//! [`WgpuRenderer`]: ../valenx_app/wgpu_renderer/struct.WgpuRenderer.html
//!
//! ## v1 strategy
//!
//! Each of the 40 functions either:
//!
//! 1. Maps cleanly onto an existing `valenx-viz` / `valenx-app` egui +
//!    wgpu capability, in which case it implements honestly (parameter
//!    validation + delegation + typed errors); or
//! 2. Returns [`OcctVizError::NotYetImplemented`] with rustdoc
//!    describing what the real OCCT API does, why direct adaptation
//!    needs more pipeline work, and which Phase `N.5` follow-up will
//!    deliver the deep implementation.
//!
//! Both kinds carry a public function signature so downstream crates
//! (the toolbox UI, the feature tree, integration tests) can be
//! written against the final API today and have the stubs fill in
//! later without churn.
//!
//! ## Feature catalogue (Phases 161-200)
//!
//! ### Viewer + camera (Phases 161-170)
//! - [`v3d_view_camera_orbit()`],
//!   [`v3d_view_camera_pan()`],
//!   [`v3d_view_camera_zoom()`],
//!   [`v3d_view_camera_perspective_toggle()`],
//!   [`v3d_view_camera_fit_all()`],
//!   [`v3d_view_camera_axo_axonometric()`],
//!   [`v3d_viewer_z_buffer()`],
//!   [`v3d_viewer_back_face_culling()`],
//!   [`v3d_viewer_clipping_plane()`],
//!   [`v3d_viewer_xor_drag()`].
//!
//! ### Selection + AIS (Phases 171-180)
//! - [`ais_interactive_context()`],
//!   [`ais_select_single()`],
//!   [`ais_select_box()`],
//!   [`ais_select_polygon()`],
//!   [`ais_highlight_dynamic()`],
//!   [`ais_highlight_selection()`],
//!   [`ais_select_face()`],
//!   [`ais_select_edge()`],
//!   [`ais_select_vertex()`],
//!   [`ais_select_subshape_filter()`].
//!
//! ### Display attributes (Phases 181-188)
//! - [`prs3d_drawer_material_default()`],
//!   [`prs3d_drawer_face_color()`],
//!   [`prs3d_drawer_edge_color()`],
//!   [`prs3d_drawer_transparency()`],
//!   [`prs3d_drawer_line_width()`],
//!   [`prs3d_drawer_line_style()`],
//!   [`prs3d_drawer_hidden_line_display()`],
//!   [`prs3d_drawer_isolines()`].
//!
//! ### Animation + transformation (Phases 189-193)
//! - [`view_animation_camera_path()`],
//!   [`view_animation_object_motion()`],
//!   [`view_animation_explode()`],
//!   [`transformation_local_axis_widget()`],
//!   [`transformation_rotation_widget()`],
//!   [`transformation_assembly_gizmo()`] — the assembly-constraint-aware
//!   gizmo: drag one part and re-solve the assembly's mates / joints so
//!   the mated parts follow (a `valenx-assembly` ↔ gizmo integration).
//!
//! ### Misc viz (Phases 194-200)
//! - [`view_screenshot()`],
//!   [`view_video_export()`],
//!   [`prs3d_drawer_section_plane_display()`],
//!   [`view_navigation_cube()`],
//!   [`view_legend_display()`],
//!   [`view_axes_origin_marker()`],
//!   [`view_grid_floor_xy()`].
//!
//! ## Error model
//!
//! All public APIs return [`Result<_, OcctVizError>`]. See
//! [`OcctVizError::code`] / [`OcctVizError::category`] for the stable
//! taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;

// Viewer + camera (Phases 161-170)
pub mod v3d_view_camera_orbit;
pub mod v3d_view_camera_pan;
pub mod v3d_view_camera_zoom;
pub mod v3d_view_camera_perspective_toggle;
pub mod v3d_view_camera_fit_all;
pub mod v3d_view_camera_axo_axonometric;
pub mod v3d_viewer_z_buffer;
pub mod v3d_viewer_back_face_culling;
pub mod v3d_viewer_clipping_plane;
pub mod v3d_viewer_xor_drag;

// Selection + AIS (Phases 171-180)
pub mod ais_interactive_context;
pub mod ais_move_to_hover;
pub mod ais_select_single;
pub mod ais_select_box;
pub mod ais_select_polygon;
pub mod ais_highlight_dynamic;
pub mod ais_highlight_selection;
pub mod ais_select_face;
pub mod ais_select_edge;
pub mod ais_select_vertex;
pub mod ais_select_subshape_filter;

// Display attributes (Phases 181-188)
pub mod prs3d_drawer_material_default;
pub mod prs3d_drawer_face_color;
pub mod prs3d_drawer_edge_color;
pub mod prs3d_drawer_transparency;
pub mod prs3d_drawer_line_width;
pub mod prs3d_drawer_line_style;
pub mod prs3d_drawer_hidden_line_display;
pub mod prs3d_drawer_isolines;

// Animation + transformation (Phases 189-193)
pub mod view_animation_camera_path;
pub mod view_animation_object_motion;
pub mod view_animation_explode;
pub mod transformation_local_axis_widget;
pub mod transformation_rotation_widget;
// Assembly-constraint-aware gizmo — drag one part, re-solve the
// assembly's mates / joints (the Tier-3 constraint-propagating drag).
pub mod transformation_assembly_gizmo;

// Misc viz (Phases 194-200)
pub mod view_screenshot;
pub mod view_video_export;
pub mod prs3d_drawer_section_plane_display;
pub mod view_navigation_cube;
pub mod view_legend_display;
pub mod view_axes_origin_marker;
pub mod view_grid_floor_xy;

pub use error::{ErrorCategory, OcctVizError};

// Re-export the entry points so callers can
// `use valenx_occt_viz::ais_select_single;` instead of the full
// module-path mouthful.
pub use ais_highlight_dynamic::ais_highlight_dynamic;
pub use ais_highlight_selection::ais_highlight_selection;
pub use ais_interactive_context::{
    ais_interactive_context, InteractiveContext, ObjectState, Pickable, PickHit, PickView,
    Ray,
};
pub use ais_move_to_hover::{move_to, HoverPreview, HoverTarget, SubKind};
pub use ais_select_box::ais_select_box;
pub use ais_select_edge::ais_select_edge;
pub use ais_select_face::ais_select_face;
pub use ais_select_polygon::ais_select_polygon;
pub use ais_select_single::ais_select_single;
pub use ais_select_subshape_filter::ais_select_subshape_filter;
pub use ais_select_vertex::ais_select_vertex;
pub use prs3d_drawer_edge_color::prs3d_drawer_edge_color;
pub use prs3d_drawer_face_color::prs3d_drawer_face_color;
pub use prs3d_drawer_hidden_line_display::prs3d_drawer_hidden_line_display;
pub use prs3d_drawer_isolines::prs3d_drawer_isolines;
pub use prs3d_drawer_line_style::prs3d_drawer_line_style;
pub use prs3d_drawer_line_width::prs3d_drawer_line_width;
pub use prs3d_drawer_material_default::prs3d_drawer_material_default;
pub use prs3d_drawer_section_plane_display::prs3d_drawer_section_plane_display;
pub use prs3d_drawer_transparency::prs3d_drawer_transparency;
pub use transformation_assembly_gizmo::{
    apply_constraint_drag, transformation_assembly_gizmo, AssemblyDragSession,
    ConstraintDragResult, DragDelta,
};
pub use transformation_local_axis_widget::{
    transformation_local_axis_widget, GizmoAxis, GizmoPlane, TranslationGizmo,
};
pub use transformation_rotation_widget::{transformation_rotation_widget, RotationGizmo};
pub use v3d_view_camera_axo_axonometric::v3d_view_camera_axo_axonometric;
pub use v3d_view_camera_fit_all::v3d_view_camera_fit_all;
pub use v3d_view_camera_orbit::v3d_view_camera_orbit;
pub use v3d_view_camera_pan::v3d_view_camera_pan;
pub use v3d_view_camera_perspective_toggle::v3d_view_camera_perspective_toggle;
pub use v3d_view_camera_zoom::v3d_view_camera_zoom;
pub use v3d_viewer_back_face_culling::v3d_viewer_back_face_culling;
pub use v3d_viewer_clipping_plane::{
    clip_mesh, v3d_viewer_clipping_plane, ClipPlane, ClipPlaneSet, MAX_CLIP_PLANES,
};
pub use v3d_viewer_xor_drag::v3d_viewer_xor_drag;
pub use v3d_viewer_z_buffer::v3d_viewer_z_buffer;
pub use view_animation_camera_path::view_animation_camera_path;
pub use view_animation_explode::view_animation_explode;
pub use view_animation_object_motion::view_animation_object_motion;
pub use view_axes_origin_marker::view_axes_origin_marker;
pub use view_grid_floor_xy::view_grid_floor_xy;
pub use view_legend_display::view_legend_display;
pub use view_navigation_cube::view_navigation_cube;
pub use view_screenshot::view_screenshot;
pub use view_video_export::{
    encode_avi, encode_h264_mp4_command, ffmpeg_available, run_ffmpeg_mp4,
    video_output_path, view_video_export, write_avi, VideoFormat, VideoFrame,
};
