use raptorq::{Decoder, Encoder, ObjectTransmissionInformation, EncodingPacket, PayloadId};
use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::{SystemTime, UNIX_EPOCH};

fn pseudo_rand_u64(state: &mut u64) -> u64 {
    // xorshift64
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage:");
        eprintln!("  encode <infile> <outfile> <packet_size> [loss] [overhead]");
        eprintln!("  decode <infile> <outfile>");
        return;
    }

    let mode = &args[1];
    let infile = &args[2];
    let outfile = &args[3];

    let mut rng_state: u64 = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(12345);
    // Mix the seed
    {
        let mut h = DefaultHasher::new();
        rng_state.hash(&mut h);
        rng_state = h.finish();
    }

    if mode == "encode" {
        let packet_size: u16 = args[4].parse().unwrap();
        let loss: f64 = if args.len() > 5 { args[5].parse().unwrap() } else { 0.0 };
        let overhead: u32 = if args.len() > 6 { args[6].parse().unwrap() } else { 5 };

        let mut data = Vec::new();
        File::open(infile).unwrap().read_to_end(&mut data).unwrap();

        let encoder = Encoder::with_defaults(&data, packet_size);
        let mut out = File::create(outfile).unwrap();
        let config = encoder.get_config();

        let f = config.transfer_length() as u64;
        let t = (config.symbol_size() - 1) as u64; 
        let oti_common = (f << 24) | (t & 0xffff);
        
        let z = (config.source_blocks() as u32) - 1;
        let n = (config.sub_blocks() as u32) - 1;
        let al = config.symbol_alignment() as u32;
        let oti_scheme = (z << 24) | (n << 8) | al;

        out.write_all(&oti_common.to_ne_bytes()).unwrap();
        out.write_all(&oti_scheme.to_ne_bytes()).unwrap();

        for block in encoder.get_block_encoders() {
            let mut num_dropped = 0u32;
            let mut kept_source_packets = Vec::new();
            
            for packet in block.source_packets() {
                let dropped: f64 = (pseudo_rand_u64(&mut rng_state) as f64) / (u64::MAX as f64) * 100.0;
                if dropped < loss {
                    num_dropped += 1;
                } else {
                    kept_source_packets.push(packet);
                }
            }
            
            let mut packets = kept_source_packets;
            packets.extend(block.repair_packets(0, num_dropped + overhead));
            
            for packet in packets {
                let sbn = packet.payload_id().source_block_number();
                let esi = packet.payload_id().encoding_symbol_id();
                
                let tag = (sbn as u32) << 24 | esi;
                out.write_all(&tag.to_ne_bytes()).unwrap();
                out.write_all(packet.data()).unwrap();
            }
        }
    } else if mode == "decode" {
        let mut in_file = File::open(infile).unwrap();
        let mut oti_common_bytes = [0u8; 8];
        in_file.read_exact(&mut oti_common_bytes).unwrap();
        let oti_common = u64::from_ne_bytes(oti_common_bytes);
        
        let mut oti_scheme_bytes = [0u8; 4];
        in_file.read_exact(&mut oti_scheme_bytes).unwrap();
        let oti_scheme = u32::from_ne_bytes(oti_scheme_bytes);

        let f = oti_common >> 24;
        let t = (oti_common & 0xffff) as u16 + 1;
        let z = ((oti_scheme >> 24) as u8) + 1;
        let n = ((oti_scheme >> 8) & 0xffff) as u16 + 1;
        let al = (oti_scheme & 0xff) as u8;

        let config = ObjectTransmissionInformation::new(f, t, z, n, al);
        let mut decoder = Decoder::new(config);

        let mut packet_data = vec![0u8; t as usize];
        let mut tag_bytes = [0u8; 4];
        let mut result = None;

        while in_file.read_exact(&mut tag_bytes).is_ok() {
            let tag = u32::from_ne_bytes(tag_bytes);
            if in_file.read_exact(&mut packet_data).is_err() {
                break;
            }
            
            let sbn = (tag >> 24) as u8;
            let esi = tag & 0x00ffffff;
            let packet = EncodingPacket::new(PayloadId::new(sbn, esi), packet_data.clone());
            let res = decoder.decode(packet);
            if res.is_some() {
                result = res;
                break;
            }
        }

        let mut out = File::create(outfile).unwrap();
        if let Some(res) = result {
            out.write_all(&res).unwrap();
        } else {
            if let Some(res) = decoder.get_result() {
                out.write_all(&res).unwrap();
            } else {
                eprintln!("Decode failed");
                std::process::exit(1);
            }
        }
    }
}

