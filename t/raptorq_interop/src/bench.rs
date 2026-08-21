use raptorq::{Decoder, Encoder, ObjectTransmissionInformation};
use std::time::Instant;
use std::env;

const TEST_BYTES: usize = 256 * 1024 * 1024; // 256 MB

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <packet_size> <num_packets> <overhead_pct>", args[0]);
        return;
    }

    let t: u16 = args[1].parse().unwrap();
    let k: usize = args[2].parse().unwrap();
    let overhead_pct: f64 = args[3].parse().unwrap();

    let expected_loss = 6.0;

    let sz = k * t as usize;
    let data: Vec<u8> = (0..sz).map(|i| (i % 256) as u8).collect();

    // ENCODE benchmark (no precalc vs precalc isn't cleanly separable in cberner API since Encoder::with_defaults does it)
    // We just measure total encode time.
    let mut bytes = 0;
    let mut elapsed_enc = 0.0;
    
    while bytes < TEST_BYTES {
        let start = Instant::now();
        let encoder = Encoder::with_defaults(&data, t);
        for block in encoder.get_block_encoders() {
            // materialize source and repair packets to simulate generation
            let _src = block.source_packets();
            let _rep = block.repair_packets(0, 0); 
        }
        elapsed_enc += start.elapsed().as_secs_f64();
        bytes += sz;
    }

    // Prepare a packet list for decode (with loss)
    let encoder = Encoder::with_defaults(&data, t);
    let mut packets_with_loss = Vec::new();
    
    for block in encoder.get_block_encoders() {
        let num_esi = block.source_packets().len();
        let mut num_dropped = 0;
        
        let mut blk_packets = Vec::new();
        for (i, p) in block.source_packets().into_iter().enumerate() {
            let dropped = ((i * 13) % 100) as f64; // simple pseudo random
            if dropped < expected_loss {
                num_dropped += 1;
            } else {
                blk_packets.push(p);
            }
        }
        let overhead = (num_esi as f64 * (overhead_pct / 100.0)) as u32;
        let rep = block.repair_packets(0, num_dropped + overhead);
        blk_packets.extend(rep);
        packets_with_loss.push(blk_packets);
    }
    
    let config = encoder.get_config();

    // DECODE benchmark (with loss)
    let mut bytes_dec = 0;
    let mut elapsed_dec = 0.0;

    while bytes_dec < TEST_BYTES {
        let start = Instant::now();
        let mut decoder = Decoder::new(config.clone());
        for blk_packets in &packets_with_loss {
            for p in blk_packets {
                decoder.decode(p.clone());
            }
        }
        let _res = decoder.get_result();
        elapsed_dec += start.elapsed().as_secs_f64();
        bytes_dec += sz;
    }

    let enc_mbps = (TEST_BYTES as f64 * 8.0) / (elapsed_enc * 1024.0 * 1024.0);
    let dec_mbps = (TEST_BYTES as f64 * 8.0) / (elapsed_dec * 1024.0 * 1024.0);

    println!("{:10} {:10.1} {:10.1}", k, enc_mbps, dec_mbps);
}

