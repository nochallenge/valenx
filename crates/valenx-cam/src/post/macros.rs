//! Macro that scaffolds a templated postprocessor.
//!
//! Each Phase 17C postprocessor variant differs from the next only in
//! its [`crate::post::template::PostTemplate`] knob values. The
//! `make_post!` macro produces a unit struct + the `Postprocessor`
//! trait impl + a simple compile-time smoke test from a template
//! description.

/// Define a unit-struct postprocessor that delegates to
/// [`crate::post::template::process_template`].
///
/// Usage:
/// ```ignore
/// make_post!(
///     name: Centroid,
///     display: "Centroid",
///     comment: CommentStyle::Open(';'),
///     prelude: "",
///     program_id: "",
///     units_mode: "G21\nG90\nG17",
///     tool_change_tpl: "T{n} M6",
///     spindle_on_tpl: "M3 S{rpm}",
///     spindle_off: "M5",
///     coolant_on: "M8",
///     coolant_off: "M9",
///     program_end: "M30",
///     program_end_suffix: "",
///     number_lines: false,
///     block_start: 10,
///     block_step: 10,
/// );
/// ```
#[macro_export]
macro_rules! make_post {
    (
        name: $struct_name:ident,
        display: $display:expr,
        comment: $comment:expr,
        prelude: $prelude:expr,
        program_id: $program_id:expr,
        units_mode: $units_mode:expr,
        tool_change_tpl: $tool_change_tpl:expr,
        spindle_on_tpl: $spindle_on_tpl:expr,
        spindle_off: $spindle_off:expr,
        coolant_on: $coolant_on:expr,
        coolant_off: $coolant_off:expr,
        program_end: $program_end:expr,
        program_end_suffix: $program_end_suffix:expr,
        number_lines: $number_lines:expr,
        block_start: $block_start:expr,
        block_step: $block_step:expr $(,)?
    ) => {
        #[doc = concat!("The ", $display, " postprocessor (Phase 17C templated variant).")]
        #[derive(Clone, Copy, Debug, Default)]
        pub struct $struct_name;

        fn tpl() -> $crate::post::template::PostTemplate {
            $crate::post::template::PostTemplate {
                name: $display,
                comment: $comment,
                prelude: $prelude,
                program_id: $program_id,
                units_mode: $units_mode,
                tool_change_tpl: $tool_change_tpl,
                spindle_on_tpl: $spindle_on_tpl,
                spindle_off: $spindle_off,
                coolant_on: $coolant_on,
                coolant_off: $coolant_off,
                program_end: $program_end,
                program_end_suffix: $program_end_suffix,
                number_lines: $number_lines,
                block_start: $block_start,
                block_step: $block_step,
            }
        }

        impl $crate::post::Postprocessor for $struct_name {
            fn header(&self, tool: &$crate::tool::Tool, spindle_rpm: f64) -> String {
                tpl().header(tool, spindle_rpm)
            }
            fn footer(&self) -> String {
                tpl().footer()
            }
            fn move_g0(&self, p: nalgebra::Vector3<f64>) -> String {
                $crate::post::format_g0(p)
            }
            fn move_g1(
                &self,
                p: nalgebra::Vector3<f64>,
                feed: f64,
            ) -> Result<String, $crate::error::CamError> {
                $crate::post::format_g1(p, feed)
            }
            fn spindle_on(&self, rpm: f64) -> String {
                tpl().spindle_on(rpm)
            }
            fn spindle_off(&self) -> String {
                tpl().spindle_off.into()
            }
            fn tool_change(&self, tool_id: u32) -> String {
                tpl().tool_change(tool_id)
            }
            fn coolant_on(&self) -> String {
                tpl().coolant_on.into()
            }
            fn coolant_off(&self) -> String {
                tpl().coolant_off.into()
            }
            fn process(
                &self,
                toolpath: &$crate::toolpath::Toolpath,
                tool: &$crate::tool::Tool,
                spindle_rpm: f64,
            ) -> Result<String, $crate::error::CamError> {
                $crate::post::template::process_template(&tpl(), toolpath, tool, spindle_rpm)
            }
        }
    };
}
