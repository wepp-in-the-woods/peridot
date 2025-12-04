use chrono::Local;

fn main() {
    // Embed a human-readable build timestamp for version reporting.
    let ts = Local::now().to_rfc3339();
    println!("cargo:rustc-env=PERIDOT_BUILD_TS={}", ts);
    println!(
        "cargo:rustc-env=PERIDOT_VERSION_STRING={} (built {})",
        env!("CARGO_PKG_VERSION"),
        ts
    );
}
