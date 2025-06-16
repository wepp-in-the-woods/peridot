use std::io::{BufRead, BufReader};
use std::fs::File;
use std::io::Write;
use std::collections::{HashSet, HashMap};

use crate::raster::Raster;

#[derive(Debug)]
pub struct ChannelNode {
    pub chnum: i32,
    pub order: i32,
    pub row: usize,
    pub col: usize,
    pub row1: usize,
    pub col1: usize,
    pub outr: usize,
    pub outc: usize,
    pub chnlen: f64,
    pub elevvup: f64,
    pub elevdn: f64,
    pub areaup: f64,
    pub areadn1: f64,
    pub areadn: f64,
    pub dda: f64,
    pub node1: i32,
    pub node2: i32,
    pub node3: i32,
    pub node4: i32,
    pub node5: i32,
    pub node6: i32,
    pub node7: i32,
    pub slopedirect: f64,
    pub slopesmoothed: f64,
    pub chn_id: i32,
    pub chnout_id: i32,
}

fn parse_line(line: &str, subwta: &Raster<i32>) -> Result<ChannelNode, Box<dyn std::error::Error>> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() != 24 {
        return Err("Incorrect number of items in line".into());
    }

    let mut node = ChannelNode {
        chnum: parts[0].parse().unwrap(),
        order: parts[1].parse().unwrap(),
        row: parts[2].parse().unwrap(),
        col: parts[3].parse().unwrap(),
        row1: parts[4].parse().unwrap(),
        col1: parts[5].parse().unwrap(),
        outr: parts[6].parse().unwrap(),
        outc: parts[7].parse().unwrap(),
        chnlen: parts[8].parse().unwrap(),
        elevvup: parts[9].parse().unwrap(),
        elevdn: parts[10].parse().unwrap(),
        areaup: parts[11].parse().unwrap(),
        areadn1: parts[12].parse().unwrap(),
        areadn: parts[13].parse().unwrap(),
        dda: parts[14].parse().unwrap(),
        node1: parts[15].parse().unwrap(),
        node2: parts[16].parse().unwrap(),
        node3: parts[17].parse().unwrap(),
        node4: parts[18].parse().unwrap(),
        node5: parts[19].parse().unwrap(),
        node6: parts[20].parse().unwrap(),
        node7: parts[21].parse().unwrap(),
        slopedirect: parts[22].parse().unwrap(),
        slopesmoothed: parts[23].parse().unwrap(),
        chn_id: -1,
        chnout_id: -1,
    };
    
    let head = (node.row - 1) * subwta.width + node.col - 1;
    let center = (node.row1 - 1) * subwta.width + node.col1 - 1;
    let tail = (node.outr - 1) * subwta.width + node.outc - 1;

    let chn_id0 = subwta.data[head];
    let chn_id1 = subwta.data[center];
    let chnout_id = subwta.data[tail];

    if chn_id0 == 0 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData,  format!("chn_id0 is zero: {}", chn_id0)))); }
    if chn_id1 == 0 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chn_id1 is zero: {}", chn_id1)))); }
    if chnout_id == 0 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chnout_id is zero: {}", chnout_id)))); }
    if chn_id0 != chn_id1 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chn_id0 != chn_id1: {} {}", chn_id0, chn_id1)))); }
    if chn_id0 % 10 != 4 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chn_id0 not ending with 4: {}", chn_id0)))); }
    if chnout_id % 10 != 4 { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chnout_id not ending with 4: {}", chnout_id)))); }
    if chn_id0 < chnout_id { return Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("chn_id0 < chnout_id: {} {}", chn_id0, chnout_id)))); }

    node.chn_id = chn_id0;
    node.chnout_id = chnout_id;

    Ok(node)
}


pub fn write_network(file_path: &str, network: &HashMap<i32, HashSet<i32>>) -> Result<(), Box<dyn std::error::Error>> {
    
    let mut network_str: Vec<String> = Vec::new();
    for (key, value) in network.iter() {
        let value_str = value.into_iter()
                             .map(|x| x.to_string())
                             .collect::<Vec<_>>()
                             .join(",");
        network_str.push(format!("{}|{}", key, value_str));
    }

    // Open a file in write mode
    let mut file = File::create(file_path)?;
    
    // Write the GeoJSON string to the file
    file.write_all(network_str.join("\n").as_bytes())?;
    Ok(())
}


pub fn read_netw_tab(file_path: &str, subwta: &Raster<i32>) -> Result<(HashMap<i32, ChannelNode>, HashMap<i32, HashSet<i32>>), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut data : HashMap<i32, ChannelNode>= HashMap::new();
    let mut network: HashMap<i32, HashSet<i32>> = HashMap::new();
    
    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        
        // Skip the first 27 lines
        if i < 27 {
            continue;
        }
        
        // Break on a line starting with "-1"
        if line.trim().starts_with("-1") {
            break;
        }
        
        // Parse line and insert data into the map
        let node = parse_line(&line, &subwta)?;

        // add chn_id0 -> chnout_id to network
        network.entry(node.chnout_id)
               .or_insert_with(HashSet::new)
               .insert(node.chn_id);
               
        data.insert(node.chn_id, node);

    }

    Ok((data, network))
}


#[cfg(test)]
mod tests {
    use std::path::Path;
    use crate::netw::read_netw_tab; 
    use crate::raster::Raster;

    #[test]
    fn test_read_netw_tab() {
        let path = "tests/fixtures/watershed_abstraction/mdobre-scarce-belch/dem/topaz/SUBWTA.ARC";
        let subwta = Raster::<i32>::read(&path).unwrap();

        let file_path = "tests/fixtures/watershed_abstraction/mdobre-scarce-belch/dem/topaz/NETW.TAB";
        match read_netw_tab(file_path, &subwta) {    
            Ok((data, network)) => {
                println!("Data: {:?}", data);
                println!("Network: {:?}", network);
        },
            Err(e) => println!("Error reading file: {}", e),
        }
    }


}
