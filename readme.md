###  Consider several possibility of the "Ns"(or rather *gaps*)
#### 1. the "Ns" can be coverd by single read, the Ns can be recovered as mismatches and be replaced, this should be most of the case for the samll gaps.
- Something has to be noted, the read has a high (0.1%, compared to the less than 0.01% mismatch and insertion) rate of deletion, so if a samll indel was found to be in a "deletion region" of the long read, the deletion should be ignored, and the gap region should be remained.
- This can be further confirmed by TRfinder, in this case, this region should be masked by tandem repeat finder, if the region is full of tandem repetor, this singnal should be igored.

#### 2. the "Ns" cause some **unaligned end**: find the begin and the end, make realignment with some aligner(bwa?muscle?), and give the Ns some repalcements, and fill the "Ns". When parsing the unaligned ends, some situations should be treated differently:
- No read support --> ignore, unfilled
- Unaligned end in both direction
  - unalign end overlapped at the end --> filled, caculated the length of the read
  - unalign end do not overlap --> double extended, unfilled
- Unaligned end in one direction
  - align with the beginning of the other edge (for example, 500bp?)
    - unalign end overlapped with the beginning of the other edge --> filled
    - unalign end do not overlap --> single extended, unfilled

#### 3. How long flaking sequence shold the be used?
- 200bp? 50bp? 20bp? For the alignment 200bp is long enough, and the pileup of samfiles can indicate the edge of the right mapping

#### 4.Keep in mind that the filling of Ns is used for the precise break point detection

