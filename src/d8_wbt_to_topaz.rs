
use crate::raster::Raster;


pub fn remap_whitebox_d8_to_topaz(flovec: &Raster<i32>) -> Raster<i32> {
    let mut remapped_flovec = flovec.empty_clone();

    for i in 0..flovec.data.len() {
        let flow_dir = flovec.data[i];
        let new_flow_dir = match flow_dir {
            1 => 3,    // East
            2 => 6,    // Northeast
            4 => 9,    // North
            8 => 8,    // Northwest
            16 => 7,   // West
            32 => 4,   // Southwest
            64 => 1,   // South
            128 => 2,  // Southeast
            _ => 0,    // No flow or undefined
        };
        remapped_flovec.data[i] = new_flow_dir;
    }

    remapped_flovec
}

// 1     2     3
//   \   |   /
//    \  | /
// 4 <-- o --> 6
//     / | \
//   /   |   \
// 7     8     9

//64	128	1
//32	0	2
//16	8	4

