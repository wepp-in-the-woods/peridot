use std::fs::OpenOptions;
use std::path::Path;
use std::sync::Once;

use env_logger::{Builder, Env, Target};

static INIT: Once = Once::new();

pub fn init_logging(log_dir: &Path) {
    INIT.call_once(|| {
        let log_path = log_dir.join("_peridot.log");
        let file = OpenOptions::new().create(true).append(true).open(&log_path);

        let mut builder = Builder::from_env(Env::default().default_filter_or("info"));
        builder.format_timestamp_millis();

        match file {
            Ok(file) => {
                builder.target(Target::Pipe(Box::new(file)));
            }
            Err(err) => {
                eprintln!("Failed to open log file {}: {}", log_path.display(), err);
                builder.target(Target::Stderr);
            }
        }

        builder.init();
    });
}
