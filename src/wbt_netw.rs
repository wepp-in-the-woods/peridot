use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
use std::fmt::Display;

#[derive(Debug)]
pub struct Link {
    pub id: i32,
    pub topaz_id: i32,
    pub ds: (isize, isize),
    pub us: (isize, isize),
    pub inflow0_id: i32,
    pub inflow1_id: i32,
    pub inflow2_id: i32,
    pub length_m: f64,
    pub ds_z: f64,
    pub us_z: f64,
    pub drop_m: f64,
    pub order: i32,
    pub areaup: f64,
    pub is_headwater: bool,
    pub is_outlet: bool,
}

fn parse_field<T>(
    s: &str,
    line_num: usize,
    name: &str
) -> io::Result<T>
where
    T: FromStr,
    <T as FromStr>::Err: Display,
{
    s.parse::<T>()
     .map_err(|e| io::Error::new(
         io::ErrorKind::InvalidData,
         format!("Failed parsing {} at line {}: {}", name, line_num+1, e)
     ))
}

pub fn read_wbt_netw_tab<P: AsRef<Path>>(file_path: P) -> 
    Result<(HashMap<i32, Link>, HashMap<i32, HashSet<i32>>), Box<dyn std::error::Error>> {


    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut links = HashMap::new();

    let mut id_to_topaz: HashMap<i32, i32> = HashMap::new();

    for (line_num, line) in reader.lines().enumerate().skip(1) { // Skip header
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        // Validate field count
        if parts.len() != 17 {
            return Err(Box::new(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Line {} has {} fields, expected 17", line_num + 1, parts.len()),
            )));
        }
        let id:        i32    = parse_field(parts[0], line_num, "id")?;
        let topaz_id:  i32    = parse_field(parts[1], line_num, "topaz_id")?;
        let ds:        (isize, isize) = {
            let x = parse_field(parts[2], line_num, "ds_x")?;
            let y = parse_field(parts[3], line_num, "ds_y")?;
            (x, y)
        };
        let us:        (isize, isize) = {
            let x = parse_field(parts[4], line_num, "us_x")?;
            let y = parse_field(parts[5], line_num, "us_y")?;
            (x, y)
        };
        let inflow0_id: i32 = parse_field(parts[6],  line_num, "inflow0_id")?;
        let inflow1_id: i32 = parse_field(parts[7],  line_num, "inflow1_id")?;
        let inflow2_id: i32 = parse_field(parts[8],  line_num, "inflow2_id")?;
        let length_m:  f64    = parse_field(parts[9],  line_num, "length_m")?;
        let ds_z:      f64    = parse_field(parts[10], line_num, "ds_z")?;
        let us_z:      f64    = parse_field(parts[11], line_num, "us_z")?;
        let drop_m:    f64    = parse_field(parts[12], line_num, "drop_m")?;
        let order:     i32    = parse_field(parts[13], line_num, "order")?;
        let areaup:    f64    = parse_field(parts[14], line_num, "areaup")?;
        let is_headwater: bool = parse_field(parts[15], line_num, "is_headwater")?;
        let is_outlet:     bool = parse_field(parts[16], line_num, "is_outlet")?;

        id_to_topaz.insert(id, topaz_id);

        // Create Link instance
        let link = Link {
            id,
            topaz_id,
            ds,
            us,
            inflow0_id,
            inflow1_id,
            inflow2_id,
            length_m,
            ds_z,
            us_z,
            drop_m,
            order,
            areaup,
            is_headwater,
            is_outlet
        };

        // Insert into HashMap
        if links.insert(topaz_id, link).is_some() {
            eprintln!("Warning: Duplicate topaz_id {} at line {}", topaz_id, line_num + 1);
        }
    }

    let mut network: HashMap<i32, HashSet<i32>> = HashMap::new();
    for link in links.values() {
        if link.inflow0_id >= 0 {
            network.entry(link.topaz_id)
                .or_insert_with(HashSet::new)
                .insert(id_to_topaz[&link.inflow0_id]);
        }
        if link.inflow1_id >= 0 {
            network.entry(link.topaz_id)
                .or_insert_with(HashSet::new)
                .insert(id_to_topaz[&link.inflow1_id]);
        }
        if link.inflow2_id >= 0 {
            network.entry(link.topaz_id)
                .or_insert_with(HashSet::new)
                .insert(id_to_topaz[&link.inflow2_id]);
        }
    }

    Ok((links, network))
}
