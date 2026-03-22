use std::collections::HashMap;
use std::fs;
use std::io::{Error, ErrorKind, Result};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct ManifestRunFlags {
    pub command: &'static str,
    pub max_points: usize,
    pub clip_hillslopes: bool,
    pub clip_hillslope_length: f64,
    pub bieger2015_widths: bool,
    pub skip_flowpaths: bool,
    pub representative_flowpath: bool,
}

#[derive(Debug, Clone)]
pub struct TabularColumn {
    pub name: &'static str,
    pub data_type: &'static str,
}

#[derive(Debug, Clone)]
pub struct TabularOutputSummary {
    pub relative_path: String,
    pub format: &'static str,
    pub rows: usize,
    pub schema: Vec<TabularColumn>,
}

fn hillslope_schema() -> Vec<TabularColumn> {
    vec![
        TabularColumn {
            name: "topaz_id",
            data_type: "int32",
        },
        TabularColumn {
            name: "slope_scalar",
            data_type: "float64",
        },
        TabularColumn {
            name: "length",
            data_type: "float64",
        },
        TabularColumn {
            name: "width",
            data_type: "float64",
        },
        TabularColumn {
            name: "direction",
            data_type: "float64",
        },
        TabularColumn {
            name: "aspect",
            data_type: "float64",
        },
        TabularColumn {
            name: "area",
            data_type: "int32",
        },
        TabularColumn {
            name: "elevation",
            data_type: "float64",
        },
        TabularColumn {
            name: "centroid_px",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_py",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_lon",
            data_type: "float64",
        },
        TabularColumn {
            name: "centroid_lat",
            data_type: "float64",
        },
    ]
}

fn channel_schema() -> Vec<TabularColumn> {
    vec![
        TabularColumn {
            name: "topaz_id",
            data_type: "int32",
        },
        TabularColumn {
            name: "slope_scalar",
            data_type: "float64",
        },
        TabularColumn {
            name: "length",
            data_type: "float64",
        },
        TabularColumn {
            name: "width",
            data_type: "float64",
        },
        TabularColumn {
            name: "direction",
            data_type: "float64",
        },
        TabularColumn {
            name: "order",
            data_type: "int32",
        },
        TabularColumn {
            name: "aspect",
            data_type: "float64",
        },
        TabularColumn {
            name: "area",
            data_type: "float64",
        },
        TabularColumn {
            name: "elevation",
            data_type: "float64",
        },
        TabularColumn {
            name: "centroid_px",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_py",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_lon",
            data_type: "float64",
        },
        TabularColumn {
            name: "centroid_lat",
            data_type: "float64",
        },
    ]
}

fn flowpath_schema() -> Vec<TabularColumn> {
    vec![
        TabularColumn {
            name: "topaz_id",
            data_type: "int32",
        },
        TabularColumn {
            name: "fp_id",
            data_type: "int32",
        },
        TabularColumn {
            name: "slope_scalar",
            data_type: "float64",
        },
        TabularColumn {
            name: "length",
            data_type: "float64",
        },
        TabularColumn {
            name: "width",
            data_type: "float64",
        },
        TabularColumn {
            name: "direction",
            data_type: "float64",
        },
        TabularColumn {
            name: "aspect",
            data_type: "float64",
        },
        TabularColumn {
            name: "area",
            data_type: "float64",
        },
        TabularColumn {
            name: "elevation",
            data_type: "float64",
        },
        TabularColumn {
            name: "order",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_px",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_py",
            data_type: "int32",
        },
        TabularColumn {
            name: "centroid_lon",
            data_type: "float64",
        },
        TabularColumn {
            name: "centroid_lat",
            data_type: "float64",
        },
    ]
}

pub fn hillslopes_summary(rows: usize, format: &'static str) -> TabularOutputSummary {
    TabularOutputSummary {
        relative_path: format!("watershed/hillslopes.{}", format),
        format,
        rows,
        schema: hillslope_schema(),
    }
}

pub fn channels_summary(rows: usize, format: &'static str) -> TabularOutputSummary {
    TabularOutputSummary {
        relative_path: format!("watershed/channels.{}", format),
        format,
        rows,
        schema: channel_schema(),
    }
}

pub fn flowpaths_summary(rows: usize, format: &'static str) -> TabularOutputSummary {
    TabularOutputSummary {
        relative_path: format!("watershed/flowpaths.{}", format),
        format,
        rows,
        schema: flowpath_schema(),
    }
}

pub fn write_watershed_readme(
    watershed_dir: &Path,
    run_flags: &ManifestRunFlags,
    tabular_outputs: &[TabularOutputSummary],
) -> Result<()> {
    if !watershed_dir.exists() {
        return Err(Error::new(
            ErrorKind::NotFound,
            format!("watershed directory missing: {}", watershed_dir.display()),
        ));
    }

    let readme_path = watershed_dir.join("README.md");
    let mut readme_size: Option<u64> = None;

    // Iterate until the README size row converges after writing.
    for _ in 0..3 {
        let content = build_readme_content(watershed_dir, run_flags, tabular_outputs, readme_size)?;
        fs::write(&readme_path, &content)?;
        let observed_size = fs::metadata(&readme_path)?.len();
        if Some(observed_size) == readme_size {
            return Ok(());
        }
        readme_size = Some(observed_size);
    }

    Ok(())
}

fn collect_files(dir: &Path, files: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            collect_files(&path, files)?;
            continue;
        }
        files.push(path);
    }
    Ok(())
}

fn to_manifest_relative(watershed_dir: &Path, file: &Path) -> Result<String> {
    let rel = file.strip_prefix(watershed_dir).map_err(|err| {
        Error::new(
            ErrorKind::InvalidInput,
            format!(
                "failed to derive relative path for {} from {}: {}",
                file.display(),
                watershed_dir.display(),
                err
            ),
        )
    })?;
    let rel_str = rel
        .to_string_lossy()
        .replace(std::path::MAIN_SEPARATOR, "/");
    Ok(format!("watershed/{}", rel_str))
}

fn detect_format(path: &Path) -> String {
    match path
        .extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or("")
        .to_ascii_lowercase()
        .as_str()
    {
        "parquet" => "parquet".to_string(),
        "csv" => "csv".to_string(),
        "geojson" => "geojson".to_string(),
        "json" => "json".to_string(),
        "txt" => "txt".to_string(),
        "slp" | "slps" => "slp".to_string(),
        "md" => "markdown".to_string(),
        ext if !ext.is_empty() => ext.to_string(),
        _ => "file".to_string(),
    }
}

fn build_readme_content(
    watershed_dir: &Path,
    run_flags: &ManifestRunFlags,
    tabular_outputs: &[TabularOutputSummary],
    readme_size: Option<u64>,
) -> Result<String> {
    let mut files: Vec<PathBuf> = Vec::new();
    collect_files(watershed_dir, &mut files)?;

    let mut manifest_rows: Vec<(String, String, u64)> = Vec::new();
    for path in files {
        let rel = to_manifest_relative(watershed_dir, &path)?;
        if rel == "watershed/README.md" {
            continue;
        }
        let format = detect_format(&path);
        let size = fs::metadata(&path)?.len();
        manifest_rows.push((rel, format, size));
    }

    manifest_rows.sort_by(|a, b| a.0.cmp(&b.0));

    let mut tabular_map: HashMap<String, &TabularOutputSummary> = HashMap::new();
    for summary in tabular_outputs {
        tabular_map.insert(summary.relative_path.clone(), summary);
    }

    let mut markdown = String::new();
    markdown.push_str("# Watershed Output Manifest\n\n");
    markdown.push_str("Generated by Peridot during watershed abstraction.\n\n");

    markdown.push_str("## Execution Flags and Inputs\n\n");
    markdown.push_str(&format!("- `command`: `{}`\n", run_flags.command));
    markdown.push_str(&format!("- `max_points`: `{}`\n", run_flags.max_points));
    markdown.push_str(&format!(
        "- `clip_hillslopes`: `{}`\n",
        run_flags.clip_hillslopes
    ));
    markdown.push_str(&format!(
        "- `clip_hillslope_length`: `{:.3}`\n",
        run_flags.clip_hillslope_length
    ));
    markdown.push_str(&format!(
        "- `bieger2015_widths`: `{}`\n",
        run_flags.bieger2015_widths
    ));
    markdown.push_str(&format!(
        "- `skip_flowpaths`: `{}`\n",
        run_flags.skip_flowpaths
    ));
    markdown.push_str(&format!(
        "- `representative_flowpath`: `{}`\n\n",
        run_flags.representative_flowpath
    ));

    markdown.push_str("## File Manifest\n\n");
    markdown.push_str("| Path | Format | Size (bytes) | Rows |\n");
    markdown.push_str("| --- | --- | ---: | ---: |\n");

    for (rel, format, size) in &manifest_rows {
        let rows = tabular_map
            .get(rel)
            .map(|summary| summary.rows.to_string())
            .unwrap_or_else(|| {
                if format == "parquet" || format == "csv" {
                    "unknown".to_string()
                } else {
                    "-".to_string()
                }
            });
        markdown.push_str(&format!("| {} | {} | {} | {} |\n", rel, format, size, rows));
    }

    let readme_rows = "-";
    let readme_size_display = readme_size
        .map(|value| value.to_string())
        .unwrap_or_else(|| "pending".to_string());
    markdown.push_str(&format!(
        "| watershed/README.md | markdown | {} | {} |\n\n",
        readme_size_display, readme_rows
    ));

    markdown.push_str("## Tabular Schema Summary\n\n");
    let mut sorted_tabular: Vec<&TabularOutputSummary> = tabular_outputs.iter().collect();
    sorted_tabular.sort_by(|a, b| a.relative_path.cmp(&b.relative_path));

    if sorted_tabular.is_empty() {
        markdown.push_str("No tabular outputs were recorded.\n\n");
    } else {
        for summary in sorted_tabular {
            markdown.push_str(&format!(
                "### `{}` (format: {}, rows: {})\n\n",
                summary.relative_path, summary.format, summary.rows
            ));
            markdown.push_str("| Column | Type |\n");
            markdown.push_str("| --- | --- |\n");
            for column in &summary.schema {
                markdown.push_str(&format!("| {} | {} |\n", column.name, column.data_type));
            }
            markdown.push('\n');
        }
    }

    markdown.push_str("## Conditional Outputs and Notes\n\n");
    if run_flags.skip_flowpaths {
        markdown.push_str(
            "- Flowpath tabular/slope outputs are skipped because `skip_flowpaths=true` (or representative mode disabled flowpaths).\n",
        );
    } else {
        markdown.push_str(
            "- Flowpath tabular and slope outputs are expected because `skip_flowpaths=false`.\n",
        );
    }

    if run_flags.representative_flowpath {
        markdown.push_str(
            "- `representative_flowpath=true` forces a single representative hillslope flowpath summary and disables full flowpath export.\n",
        );
    }

    if run_flags.clip_hillslopes {
        markdown.push_str(&format!(
            "- Hillslope slope profiles are clipped to `clip_hillslope_length={:.3}`.\n",
            run_flags.clip_hillslope_length
        ));
    } else {
        markdown
            .push_str("- Hillslope slope profiles are not clipped (`clip_hillslopes=false`).\n");
    }

    if run_flags.bieger2015_widths {
        markdown.push_str("- Channel widths use the Bieger 2015 width model.\n");
    } else {
        markdown.push_str("- Channel widths use the legacy width model.\n");
    }

    Ok(markdown)
}
