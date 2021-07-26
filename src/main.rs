use std::sync::Arc;
use std::thread;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;

use bio::io::fasta;
use bio::io::fasta::FastaRead;


fn main() {
	let mut reader = fasta::Reader::from_file("CP000648.1.fasta").unwrap();
	let mut record = fasta::Record::new();

	while let Ok(()) = reader.read(&mut record) {
		// Create an FM-Index for a given text.
		let alphabet = alphabets::dna::iupac_alphabet();

		//>NG_049244.1 Klebsiella pneumoniae Kp10-26 pKp10-26 blaKPC gene for carbapenem-hydrolyzing class A beta-lactamase KPC-11, complete CDS
		let amrs = vec![
		b"ATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCG",
		b"CGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTA",
		b"CGCGATGGATACCGGCTCAGGCGCAACTGTAAGTTACCGCGCTGAGGAGCGCTTCCCACTGTGCAGCTCA",
		b"TTCAAGGGCTTTCTTGCTGCCGCTGTGCTGGCTCGCAGCCAGCAGCAGGCCGGCTTGCTGGACACACCCA",
		b"TCCGTTACGGCAAAAATGCGCTGGTTCTGTGGTCACCCATCTCGGAAAAATATCTGACAACAGGCATGAC",
		b"GGTGGCGGAGCTGTCCGCGGCCGCCGTGCAATACAGTGATAACGCCGCCGCCAATTTGTTGCTGAAGGAG",
		b"TTGGGCGGCCCGGCCGGGCTGACGGCCTTCATGCGCTCTATCGGCGATACCACGTTCCGTCTGGACCGCT",
		b"GGGAGCTGGAGCTGAACTCCGCCATCCCAGGCGATGCGCGCGATACCTCATCGCCGCGCGCCGTGACGGA",
		b"AAGCTTACAAAAACTGACACTGGGCTCTGCACTGGCTGCGCCGCAGCGGCAGCAGTTTGTTGATTGGCTA",
		b"AAGGGAAACACGACCGGCAACCACCGCATCCGCGCGGCGGTGCCGGCAGACTGGGCAGTCGGAGACAAAA",
		b"CCGGAACCTGCGGAGTGTATGGCACGGCAAATGACTATGCCGTCGTCTGGCCCACTGGGCGCGCACCTAT",
		b"TGTGTTGGCCGTCTACACCCGGGCGCCTAACAAGGATGACAAGCACAGCGAGGCCGTCATCGCCGCTGCG"];

	    if record.is_empty() {
	        break;
	    }

		let check = record.check();
	    if check.is_err() {
	        panic!("I got a rubbish record!")
	    }
	    // obtain sequence
	    let text = record.seq();
	    // check, whether seq is in the expected alphabet
	    if alphabet.is_word(text) {
			println!("{}", record.id());
			let sa = suffix_array(text);
			let bwt = Arc::new(bwt(text, &sa));
			let less = Arc::new(less(bwt.as_ref(), &alphabet));
			let occ = Arc::new(Occ::new(bwt.as_ref(), 3, &alphabet));
			let fmindex = Arc::new(FMIndex::new(bwt, less, occ));

			// Spawn threads to perform backward searches for each interval
			let interval_calculators = amrs
			    .into_iter()
			    .map(|amr| {
			        let fmindex = fmindex.clone();
			        thread::spawn(move || fmindex.backward_search(amr.iter()))
			    })
			    .collect::<Vec<_>>();

			// Loop through the results, extracting the positions array for each pattern
			for interval_calculator in interval_calculators {
				let result = interval_calculator.join().unwrap();
			    let positions = result.occ(&sa);
				println!("{:?}", positions);
				//println!("{:?}", result);
				println!("{:?}", &text[0..4]);
			}
	    }
	    reader.read(&mut record).expect("Failed to parse record");
	}
}
