
use crate::raster::Raster;

pub fn remap_whitebox_d8_to_topaz_in_place(flovec: &mut Raster<u8>) {
    for value in flovec.data.iter_mut() {
        *value = match *value {
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
    }
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
