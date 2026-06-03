//! Persistence helpers shared between the STEP and IGES backends.
//!
//! Today this module just centralises the timestamp + filename header
//! formatting used in both IGES Global sections and STEP `FILE_NAME`
//! attributes so the two backends emit consistent metadata.

use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

/// IGES-formatted timestamp (`YYYYMMDD.HHMMSS`).
///
/// IGES 5.3 Global section parameters 18 and 19 want a 15-character
/// date/time stamp. We don't link to the `chrono` crate just for this
/// — synthesize directly from `SystemTime`.
pub fn iges_timestamp(now: SystemTime) -> String {
    let secs = now
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (year, month, day, hour, minute, second) = unix_to_ymd_hms(secs);
    format!("{year:04}{month:02}{day:02}.{hour:02}{minute:02}{second:02}")
}

/// Strip the directory portion of a path; fall back to a constant when
/// the path has no file component (e.g. `"/"`).
pub fn basename(path: &Path) -> String {
    path.file_name()
        .and_then(|s| s.to_str())
        .map(|s| s.to_owned())
        .unwrap_or_else(|| "valenx_export".to_string())
}

/// Convert a Unix timestamp (seconds since 1970-01-01 UTC) into
/// `(year, month, day, hour, minute, second)`. Pure arithmetic — no
/// timezone handling, no leap-second drift; good enough for a stamp on
/// the header of an export file.
fn unix_to_ymd_hms(secs: i64) -> (i32, u32, u32, u32, u32, u32) {
    let mut days = secs / 86_400;
    let mut rem = secs.rem_euclid(86_400) as u32;
    if secs < 0 && rem != 0 {
        days -= 1;
    }
    let hour = rem / 3_600;
    rem %= 3_600;
    let minute = rem / 60;
    let second = rem % 60;

    // Civil-from-days algorithm by Howard Hinnant (CC BY 4.0).
    let z = days + 719_468;
    let era = if z >= 0 {
        z / 146_097
    } else {
        (z - 146_096) / 146_097
    };
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let year = (y + if month <= 2 { 1 } else { 0 }) as i32;
    (year, month, day, hour, minute, second)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iges_timestamp_known_epoch() {
        // 2000-01-01 00:00:00 UTC = 946_684_800
        let t = UNIX_EPOCH + std::time::Duration::from_secs(946_684_800);
        let stamp = iges_timestamp(t);
        assert_eq!(stamp, "20000101.000000", "got {stamp}");
    }

    #[test]
    fn iges_timestamp_2026_05_17() {
        // 2026-05-17 12:34:56 UTC = 1_779_021_296
        let t = UNIX_EPOCH + std::time::Duration::from_secs(1_779_021_296);
        let stamp = iges_timestamp(t);
        assert_eq!(stamp, "20260517.123456", "got {stamp}");
    }

    #[test]
    fn basename_strips_directories() {
        assert_eq!(basename(Path::new("/tmp/foo.step")), "foo.step");
        assert_eq!(basename(Path::new("foo.iges")), "foo.iges");
    }

    #[test]
    fn basename_falls_back_on_rootless_path() {
        // A path that is just `/` has no file name — fall back to the
        // constant so the header is never blank.
        assert_eq!(basename(Path::new("/")), "valenx_export");
    }
}
