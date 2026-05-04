//! n2core/src/sequence.rs


pub trait DnaSequence {
    /// Owned type (String or Vec<u8>)
    type OwnedSeq;
    
    /// Borrowed slice (&str or &[u8])
    type SeqSlice<'a> where Self: 'a;

    /// Returns the reverse complement of the sequence
    fn reverse_complement(&self) -> Self::OwnedSeq;
    
    /// Returns an iterator yielding consecutive k-mers of length `k`
    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>>;
}

// ---------------------------------------------------------
// Implementation for str
// ---------------------------------------------------------
impl DnaSequence for str {
    type OwnedSeq = String;
    type SeqSlice<'a> = &'a str;

    fn reverse_complement(&self) -> Self::OwnedSeq {
        self.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T', 'a' => 't',
                'C' => 'G', 'c' => 'g',
                'G' => 'C', 'g' => 'c',
                'T' => 'A', 't' => 'a',
                'N' => 'N', 'n' => 'n',
                _ => c, 
            })
            .collect()
    }

    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>> {
        self.as_bytes()
            .windows(k)
            .map(|w| std::str::from_utf8(w).expect("Invalid UTF-8 in DNA string"))
    }
}

// ---------------------------------------------------------
// Implementation for u8
// ---------------------------------------------------------
impl DnaSequence for [u8] {
    type OwnedSeq = Vec<u8>;
    type SeqSlice<'a> = &'a [u8];

    fn reverse_complement(&self) -> Self::OwnedSeq {
        self.iter()
            .rev()
            .map(|&c| match c {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'N' | b'n' => b'N',
                _ => c, // Keep other bytes as-is
            })
            .collect()
    }

    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>> {
        self.windows(k)
    }
}
