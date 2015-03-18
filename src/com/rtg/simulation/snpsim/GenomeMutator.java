/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.simulation.snpsim;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.simulation.snpsim.Mutation.MutationType;
import com.rtg.util.Constants;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.VariantStatistics;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.VcfHeader;

/**
 * Mutates a genome
 */
public class GenomeMutator {
  GenomeMutatorPriors mMutationGenerator;

  // Average failed positions are you prepared to accept in order to create a new mutation
  static final int VALID_RATIO = 10;
  // Burst limit for failed mutation positions.
  static final int INVALID_THRESHOLD = 200;

  private static final byte[] MAPPING;
  static {
    MAPPING = new byte[DNA.values().length];
    for (int i = 0; i < MAPPING.length; i++) {
      MAPPING[i] = (byte) DNA.values()[i].toString().charAt(0);
    }
  }

  protected final boolean mVerbose;
  protected final PortableRandom mRandom;
  private final boolean mSimpleMnps;
  private final int mMinMutationDistance;

  private GenomePriorParams mPriors;
  protected VcfWriter mMapping;
  private VcfHeader mHeader;

  final VariantOutputVcfFormatter mFormatter;
  final VariantParams mVariantParams;
  VariantStatistics mVariantStatistics;
  final String mSample;

  // input parameter rates
  private Double mSnpRate, mMnpRate, mIndelRate;

  GenomeMutator(int seed, boolean verbose, boolean simpleMnps, int minDistance, VariantParams variantParams) {
    this(new PortableRandom(seed), verbose, simpleMnps, minDistance, variantParams, "SAMPLE");
  }


  GenomeMutator(PortableRandom random, boolean verbose, boolean simpleMnps, int minDistance, VariantParams variantParams, String sample) {
    mSample = sample;
    mRandom = random;
    mVerbose = verbose;
    mSimpleMnps = simpleMnps;
    mMinMutationDistance = minDistance;
    mFormatter = new VariantOutputVcfFormatter(variantParams, sample);
    mVariantStatistics = new VariantStatistics(variantParams.directory());
    mVariantStatistics.onlySamples(mSample);
    mVariantParams = variantParams;
  }

  protected static MutationComparator getMutationComparator() {
    return new MutationComparator();
  }

  protected static class MutationComparator implements Comparator<Mutation>, Serializable {
    @Override
    public int compare(final Mutation o1, final Mutation o2) {
      return Integer.valueOf(o1.mPos).compareTo(o2.mPos);
    }
  }

  void setPriors(final GenomePriorParams pp) {
    mPriors = pp;
  }
  void setRates(final Double snpRate, final Double mnpRate, final Double indelRate) {
    mSnpRate = snpRate;
    mMnpRate = mnpRate;
    mIndelRate = indelRate;
  }

  void setSimplePriors(final double mRate) {
    mMutationGenerator = new GenomeMutatorPriors(mSnpRate, mMnpRate, mIndelRate, mRate);
  }

  int mutatePriors(final SequencesReader input, Sex sex, final File outputDirectory, final File twinDirectory, final OutputStream mappingOutput) throws IOException {
    assert mPriors != null;
    mMutationGenerator = new GenomeMutatorPriors(mPriors, twinDirectory == null ? Ploidy.HAPLOID : Ploidy.DIPLOID);
    final double rate = mMutationGenerator.rate();
    return mutateRate(input, sex, outputDirectory, twinDirectory, mappingOutput, rate);
  }

  int mutateCount(final SequencesReader input, Sex sex, final File outputDirectory, final File twinDirectory, final OutputStream mappingOutput, final int count) throws IOException {
    assert mPriors == null;
    final double rate = (double) count / input.totalLength();
    if (rate > 1.0) {
      Diagnostic.warning("Count gives a mutation rate of greater than 1 mutation per nucleotide");
      return 1;
    }
    return mutateRate(input, sex, outputDirectory, twinDirectory, mappingOutput, rate);
  }

  int mutateRate(SequencesReader input, Sex sex, File outputDirectory, File twinDirectory, OutputStream mappingOutput, double rate) throws IOException {
    if (rate * input.maxLength() < 1) {
      Diagnostic.warning("Mutation rate is less than 1 mutation per sequence");
    }
    if (mPriors == null) {
      setSimplePriors(rate);
    }
    final PrereadType pType = input.getPrereadType();
    final ReferenceGenome rg = new ReferenceGenome(input, sex, DefaultFallback.DIPLOID); //TODO: allow for HAPLOID Default also
    final SdfWriter output = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE, pType, false, true, true, input.type());
    try {
      output.setPrereadArm(input.getArm());
      output.setCommandLine(CommandLine.getCommandLine());

      final SdfWriter twinOutput;
      if (twinDirectory != null) {
        twinOutput = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE, pType, false, true, true, input.type());
        twinOutput.setPrereadArm(input.getArm());
        twinOutput.setCommandLine(CommandLine.getCommandLine());
        twinOutput.setComment("Sex=" + sex + " HaploidSet=FIRST");
        output.setComment("Sex=" + sex + " HaploidSet=SECOND");
      } else {
        twinOutput = null;
      }
      try {
        mHeader = mFormatter.makeHeader(mVariantParams, null);
        mMapping = new VcfWriter(mHeader, mappingOutput);
        try {
          for (int i = 0; i < input.numberSequences(); i++) {
            final byte[] genomeSeq = new byte[input.length(i)];
            final String seqName = input.name(i);
            final ReferenceSequence refSeq = rg.sequence(seqName);
            final Ploidy ploidy = twinOutput == null ? Ploidy.HAPLOID : refSeq.ploidy();
            if (ploidy == Ploidy.NONE) { // This sequence is not present in this sex of individual
              continue;
            }
            final int length = input.read(i, genomeSeq);
            final Mutation[] sortedLocs = generatePositions(rate, ploidy, genomeSeq, length);
            final SdfWriter h1;
            final SdfWriter h2;
            if (ploidy == Ploidy.HAPLOID && twinOutput != null && refSeq.haploidComplementName() != null) {
              if (refSeq.hasDuplicates()) {
                Diagnostic.warning("PAR regions not supported, ignoring.");
              }
              final String hapName = refSeq.haploidComplementName();
              // For a haploid complement, work out which destination gets the sequence
              if (hapName.compareTo(seqName) > 0) {
                h1 = output;
              } else {
                h1 = twinOutput;
              }
              h2 = null;
            } else {
              h1 = output;
              h2 = twinOutput;
            }
            h1.startSequence(seqName);
            if (h2 != null) {
              h2.startSequence(seqName);
            }
            mutateSequence(h1, h2, genomeSeq, sortedLocs, seqName);
            h1.endSequence();
            if (h2 != null) {
              h2.endSequence();
            }
          }
        } finally {
          mMapping.close();
        }
      } finally {
        if (twinOutput != null) {
          twinOutput.close();
        }
      }
    } finally {
      output.close();
    }
    return 0;
  }

  void summarise(File outputDir, File twinOutputDir, OutputStream out) throws IOException {
    final byte[] stats = statisticsLines().getBytes();
    final File summary = new File(outputDir, CommonFlags.SUMMARY_FILE);
    try (FileOutputStream fos = new FileOutputStream(summary)) {
      fos.write(stats);
    }
    if (twinOutputDir != null) {
      final File twinsummary = new File(twinOutputDir, CommonFlags.SUMMARY_FILE);
      try (FileOutputStream tfos = new FileOutputStream(twinsummary)) {
        tfos.write(stats);
      }
    }
    out.write(stats);
  }

  /**
   * @param rate rate of mutations
   * @param ploidy the ploidy of this genome sequence
   * @param genomeSeq template sequence
   * @param length the length of the template sequence
   * @return a list of mutations
   * @throws IOException
   */
  private Mutation[] generatePositions(final double rate, Ploidy ploidy, final byte[] genomeSeq, int length) throws IOException {
    final int mutationNo = (int) (rate * length);
    final Mutation[] locs = new Mutation[mutationNo];
    int i = 0;
    int trying = 0;
    // will terminate if there are more than 10 times as many failures
    // as successes for an extended period
    while (i < mutationNo && trying < INVALID_THRESHOLD) {
      final int pos = mRandom.nextInt(length);
      final Mutation mut = new Mutation(mRandom, pos, mMutationGenerator, length, genomeSeq[pos], ploidy);
      if (validPos(locs, i, mut, genomeSeq)) {
        locs[i] = mut;
        if (ploidy == Ploidy.HAPLOID) {
          locs[i] = reduceToHomozygous(locs[i]);
        }
        trying = Math.max(0, trying - VALID_RATIO);
        i++;
      } else {
        trying++;
      }
    }
    final Mutation[] sortedLocs = new Mutation[i];
    System.arraycopy(locs, 0, sortedLocs, 0, i);
    if (sortedLocs.length > 0) {
      Arrays.sort(sortedLocs, getMutationComparator());
    }
    return sortedLocs;
  }

  private Mutation reduceToHomozygous(final Mutation m) {
    Mutation ret = m;
    if (m.mHeterozygous ||  m.mGenDiffMode != GenDiffMode.BOTH_SAME) {
      ret = new Mutation(m.mPos, m.mType, false, m.mDiffMode, GenDiffMode.BOTH_SAME,
          m.mLength != 0 ? m.mLength : m.mLengthTwin, m.mLength != 0 ? m.mLength : m.mLengthTwin, m.mBases, Ploidy.HAPLOID);
    }
    return ret;
  }

  // if not fast enough will have to be improved  sofar does not allow SNPs into diploid-deletes
  protected boolean validPos(final Mutation[] locs, final int currentNumLocs, Mutation mut, byte[] seq) {
    final int pos = mut.mPos;
    final int maxMutation;
    if (mut.mType != MutationType.INSERT) {
      maxMutation = mut.mLength > mut.mLengthTwin ? mut.mLength : mut.mLengthTwin;
      // check for Ns
      int i = 0;
      while (i < maxMutation) {
        if (!getNucleotideString(seq, pos + i, 1).equals("N")) {
          break;
        }
        i++;
      }
      if (i == maxMutation) {
        return false;
      }
    } else {
      if (getNucleotideString(seq, pos, 1).equals("N")) {
        return false;
      }
      maxMutation = 0;
    }
    // check for position overlap
    for (int j = 0; j < currentNumLocs; j++) {
      final int maxSeq = locs[j].mLength > locs[j].mLengthTwin ? locs[j].mLength : locs[j].mLengthTwin;
      switch (locs[j].mType.simple()) {
      case SNP:
      case MNP:
      case DELETE:
        if (pos + maxMutation - 1 + mMinMutationDistance >= (locs[j].mPos)
            && pos < (locs[j].mPos + maxSeq  - 1 + mMinMutationDistance)) {
          return false;
        }
        break;
      case INSERT:
      case INSDEL:
        if (pos > (locs[j].mPos - mMinMutationDistance)
            && pos < (locs[j].mPos + maxSeq + mMinMutationDistance)) {
          return false;
        }
        break;
      default:
        throw new IllegalStateException("Unpossible");
      }
    }
    return true;
  }

  int mutateMnp(SdfWriter h1, SdfWriter h2, final byte[] seq, final String seqName, final Mutation m) throws IOException {
    final int lenMax = m.mLength > m.mLengthTwin ? m.mLength : m.mLengthTwin;
    final String mut1;
    final String mut2;
    final String old = getNucleotideString(seq, m.mPos, lenMax);

    DNA[] dnaArray1 = new DNA[0]; DNA[] dnaArray2 = new DNA[0];
    if (m.mType == MutationType.SNP && mPriors != null) {
      mut1 = "" + DnaUtils.getBase(m.mBases[0]);
      mut2 = "" + DnaUtils.getBase(m.mBases[1]);
    } else {
      if (m.mLength > 0) {
        if (mSimpleMnps) {
          dnaArray1 = getMutatedMnpSimple(seq, getDnaArray(seq, m.mPos, m.mLength), m.mPos, m.mLength, m.mLengthTwin);
        } else {
          dnaArray1 = getMutatedMnp(seq, getDnaArray(seq, m.mPos, m.mLength), m.mPos, m.mLength, m.mLengthTwin);
        }
      }
      if (m.mGenDiffMode == GenDiffMode.TWIN_ONLY && m.mLengthTwin > 0) {
        if (mSimpleMnps) {
          dnaArray2 = getMutatedMnpSimple(seq, getDnaArray(seq, m.mPos, m.mLengthTwin), m.mPos, m.mLengthTwin, m.mLength);
        } else {
          dnaArray2 = getMutatedMnp(seq, getDnaArray(seq, m.mPos, m.mLengthTwin), m.mPos, m.mLengthTwin, m.mLength);
        }
      }
      if (m.mLengthTwin > 0 && m.mGenDiffMode == GenDiffMode.DIFFERENT) {
        if (mSimpleMnps) {
          dnaArray2 = getMutatedMnpSimple(seq, dnaArray1, m.mPos, m.mLength, m.mLengthTwin);
        } else {
          dnaArray2 = getMutatedMnp(seq, dnaArray1, m.mPos, m.mLengthTwin, m.mLength);
        }
      }
      mut1 = dnaArray1.length == 0 ? "i" : dnasToString(dnaArray1);
      mut2 = dnaArray2.length == 0 ? "i" : dnasToString(dnaArray2);
    }

    final VcfRecord rec;
    if (m.mHeterozygous) {
      final String[] ploids = GenDiffMode.diffModeArray(new String[] {old, mut1, mut2}, m.mGenDiffMode);
      rec = mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, ploids);
    } else { // homozygous
      rec = mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, old, mut1);
    }
    if (mMapping != null) { // Stupid tests do this
      mMapping.write(rec);
    }
    return lenMax;
  }

  protected DNA[] getMutatedMnp(final byte[] seq, final DNA[] oldDna2, final int pos, final int length, final int otherLength) {
    final boolean snpOnly = length <= 2 && otherLength <= 2;
    final int maxLength = length > otherLength ? length : otherLength;

    DNA[] dnaArray;
    do {
      final ArrayList<DNA[]> dnaList = new ArrayList<>(maxLength);
      int newLength = 0;
      //                                                SNP     ins  delete NOP
      final double[] midThresholds = nextBoolean(mRandom) ? GenomeMutatorPriors.THRESHOLD_MID1 : GenomeMutatorPriors.THRESHOLD_MID2;
      final double[] endThresholds = snpOnly ? GenomeMutatorPriors.THRESHOLD_ENDS_SNPS : GenomeMutatorPriors.THRESHOLD_ENDS;
      for (int i = 0; i < length; i++) {
        final double[] thresholds = (i == 0 || i == length - 1) ? endThresholds : midThresholds;
        final DNA[] newDna = chooseMnpMutation(mRandom.nextDouble(), DNA.values()[seq[pos + i]], oldDna2.length > i ? oldDna2[i] : DNA.N, thresholds);
        dnaList.add(newDna);
        newLength += newDna.length;
      }
      final int rest = maxLength - length;
      dnaArray = new DNA[newLength + rest];
      int index = 0;
      for (final DNA[] dnaPiece : dnaList) {
        for (final DNA x : dnaPiece) {
          dnaArray[index++] = x;
        }
      }
      for (int i = length; i < maxLength; i++) {
        dnaArray[index++] = DNA.values()[seq[pos + i]];
      }
    } while (dnasToString(dnaArray).equals(dnasToString(oldDna2)));
    return dnaArray;
  }

  /* Generates a mutated MNP of same length with each position a SNP */
  private DNA[] getMutatedMnpSimple(final byte[] seq, final DNA[] oldDna2, final int pos, final int len, final int totalLen) {
    final int lenMax = len > totalLen ? len : totalLen;
    final DNA[] dnaArray = new DNA[lenMax];
    for (int i = 0; i < len; i++) {
      final DNA old = DNA.values()[seq[pos + i]];
      final DNA old2 = oldDna2[i];
      DNA mut1 = randomDNA(mRandom);
      while (mut1 == old || mut1 == old2) {
        mut1 = randomDNA(mRandom);
      }
      dnaArray[i] = mut1;
    }
    for (int i = len; i < totalLen; i++) {
      dnaArray[i] = DNA.values()[seq[pos + i]];
    }
    return dnaArray;
  }

  DNA[] chooseMnpMutation(double rand, DNA old, DNA old2, double[] thres) {
    final DNA[] newDna;
    if (rand < thres[0]) {
      // snp
      DNA mut1 = randomDNA(mRandom);
      while (mut1 == old || mut1 == old2) {
        mut1 = randomDNA(mRandom);
      }
      newDna = new DNA[1];
      newDna[0] = mut1;
    } else if (rand < thres[1]) {
      // insert
      final int len = mMutationGenerator.chooseLength(mRandom.nextDouble(), MutationType.INSERT, true);
      newDna = new DNA[len + 1];
      for (int i = 0; i < newDna.length - 1; i++) {
        newDna[i] = randomDNA(mRandom);
      }
      newDna[newDna.length - 1] = old;
    } else if (rand < thres[2]) {
      // delete
      newDna = new DNA[0];
    } else {
      newDna = new DNA[1];
      newDna[0] = old;
    }
    return newDna;
  }

  static String dnasToString(final DNA[] dnas) {
    if (dnas == null) {
      throw new IllegalArgumentException("MNP dnas shouldnt be null");
    }
    final byte[] characters = new byte[dnas.length];
    for (int i = 0; i < dnas.length; i++) {
      characters[i] = (byte) dnas[i].toString().charAt(0);
    }
    return new String(characters);
  }

  private int mutateDelete(SdfWriter h1, SdfWriter h2, final byte[] seq, final String seqName, final Mutation m) throws IOException {
    assert m.mType == Mutation.MutationType.DELETE;
    final int maxlen = m.mLength > m.mLengthTwin ? m.mLength : m.mLengthTwin;
    final String old = getNucleotideString(seq, m.mPos, maxlen);
    assert !(h2 == null && m.mLengthTwin != m.mLength);
    final String[] mutations = generateDeletion(m, old, mRandom);
    assert !(h2 == null && m.mGenDiffMode != GenDiffMode.BOTH_SAME);
    assert !(m.mGenDiffMode != GenDiffMode.BOTH_SAME && m.mLength == m.mLengthTwin);
    final VcfRecord rec = mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, mutations);
    if (mMapping != null) { // Stupid unit tests do this
      mMapping.write(rec);
    }
    return maxlen;
  }

  /**
   * @param m the mutation specification
   * @param old the original nucleotides in this region
   * @param random a source of randomness
   * @return an array containing the mutation triplet (old, new first, new second)
   */
  protected static String[] generateDeletion(final Mutation m, final String old, final PortableRandom random) {
    final int difflen = Math.abs(m.mLength - m.mLengthTwin);
    assert old.length() == Math.max(m.mLength, m.mLengthTwin);
    final StringBuilder ds1 = new StringBuilder();
    final StringBuilder ds2 = new StringBuilder();
    if (m.mLength < m.mLengthTwin) {
      ds2.append('i');
      if (m.mLength == 0) {
        ds1.append(old);
      } else {
        final int delPos = random.nextInt(difflen);
        ds1.append(old.substring(0, delPos));
        ds1.append(old.substring(delPos + m.mLength, old.length()));
      }
    } else if (m.mLengthTwin < m.mLength) {
      ds1.append('i');
      if (m.mLengthTwin == 0) {
        ds2.append(old);
      } else {
        final int delPos = random.nextInt(difflen);
        ds2.append(old.substring(0, delPos));
        ds2.append(old.substring(delPos + m.mLengthTwin, old.length()));
      }
    } else {
      ds1.append('i');
      ds2.append('i');
    }
    return GenDiffMode.diffModeArray(new String[] {old, ds1.toString(), ds2.toString()}, m.mGenDiffMode);
  }

  static DNA randomDNA(final PortableRandom r) {
    return DNA.values()[1 + r.nextInt(4)];
  }

  private int mutateInsert(SdfWriter h1, SdfWriter h2, final byte[] seq, final String seqName, final Mutation m) throws IOException {
    final String[] mutationParts = generateInsertion(m, mRandom);
    final VcfRecord rec = mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, mutationParts);
    if (mMapping != null) { // Stupid unit tests may do this
      mMapping.write(rec);
    }
    return 0;
  }

  /**
   * @param m the mutation specifying position in genome to generate
   * @param random a source of randomness
   * @return an array containing up to three strings specifying the insert.
   */
  static String[] generateInsertion(final Mutation m, final PortableRandom random) {
    assert m.mType == MutationType.INSERT;
    //two insertions
    final int len1 = m.mLength > 0 ? m.mLength : m.mLengthTwin;
    final DNA[] dnaArray1 = new DNA[len1];
    for (int i = 0; i < len1; i++) {
      dnaArray1[i] = randomDNA(random);
    }
    final String[] mutationParts;
    if (m.mHeterozygous) {
      // need a second dna string
      final int len2 = m.mLengthTwin;
      final DNA[] dnaArray2 = new DNA[len2];
      boolean done = false;
      while (!done) {
        for (int i = 0; i < len2; i++) {
          dnaArray2[i] = randomDNA(random);
        }
        if (len1 != len2) {
          done = true;
        } else {
          for (int i = 0; i < len2 ; i++) {
            if (dnaArray2[i] != dnaArray1[i]) {
              done = true;
              break;
            }
          }
        }
      }
      final String[] parts = {"i", dnasToString(dnaArray1), dnasToString(dnaArray2)};
      mutationParts = GenDiffMode.diffModeArray(parts, m.mGenDiffMode);
    } else { // homozygous case
      mutationParts = new String[] {"i", dnasToString(dnaArray1)};
    }
    return mutationParts;
  }

  private int mutateInsDel(SdfWriter h1, SdfWriter h2, final byte[] seq, final String seqName, final Mutation m) throws IOException {
    assert m.mType == MutationType.INSDEL;
    // make insert and deletion, decide (50/50) which one first
    final boolean delFirst = nextBoolean(mRandom);
    // prepare deletion
    final int delLen = delFirst ? m.mLength : m.mLengthTwin;
    final String old = getNucleotideString(seq, m.mPos, delLen);

    final int insLen = delFirst ? m.mLengthTwin : m.mLength;
    final DNA[] dnaArray = new DNA[insLen];
    for (int i = 0; i < insLen; i++) {
      dnaArray[i] = randomDNA(mRandom);
    }
    if (m.mHeterozygous) {
      if (delFirst) {
        mMapping.write(mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos),  mVerbose, old, "i", dnasToString(dnaArray) + old));
      } else {
        mMapping.write(mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, old, dnasToString(dnaArray) + old, "i"));
      }
    } else { // homozygous case
      if (delFirst) {
        mMapping.write(mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, old, "i"));
      } else {
        mMapping.write(mutationLine(h1, h2, seqName, m, getPrevRefChar(seq, m.mPos), mVerbose, old, dnasToString(dnaArray)));
      }
    }
    return delLen;  // jump is delete length
  }

  private char getPrevRefChar(byte[] ref, int pos) {
    return pos <= 0 ? 'N' : DnaUtils.getBase(ref[pos - 1]);
  }

  private boolean nextBoolean(PortableRandom random) {
    return random.nextDouble() > 0.5;
  }

  protected void mutateSequence(SdfWriter h1, SdfWriter h2, final byte[] seq, final Mutation[] mutations, final String seqName) throws IOException {
    int i = 0;
    int pos = 0;
    while (pos < seq.length) {
      //write out regular tides until location
      final int end = i < mutations.length ? mutations[i].mPos : seq.length;
      if (end - pos > 0) {
        final byte[] buf = new byte[end - pos];
        System.arraycopy(seq, pos, buf, 0, buf.length);
        if (h1 != null) {
          h1.write(buf, null, buf.length);
        }
        if (h2 != null) {
          h2.write(buf, null, buf.length);
        }
        pos = end;
      }
      if (i < mutations.length) {
        if (mutations[i].mLength > 0 || mutations[i].mLengthTwin > 0) {
          switch (mutations[i].mType.simple()) {
          case SNP:
          case MNP:
            pos += mutateMnp(h1, h2, seq, seqName, mutations[i]);
            break;
          case INSERT:
            pos += mutateInsert(h1, h2, seq, seqName, mutations[i]);
            break;
          case DELETE:
            pos += mutateDelete(h1, h2, seq, seqName, mutations[i]);
            break;
          case INSDEL:
            pos += mutateInsDel(h1, h2, seq, seqName, mutations[i]);
            break;
          default:
            throw new IllegalStateException("Unpossible");
          }
        }
        i++;
      }
    }
  }

  static String getNucleotideString(final byte[] seq, final int pos, final int length) {
    final byte[] ret = new byte[length];
    for (int i = 0; i < length; i++) {
      ret[i] = MAPPING[seq[pos + i]];
    }
    return new String(ret);
  }

  private DNA[] getDnaArray(final byte[] seq, final int pos, final int length) {
    final DNA[] arr = new DNA[length];
    for (int i = 0; i < length; i++) {
      arr[i] = DNA.values()[seq[pos + i]];
    }
    return arr;
  }

  String statisticsLines() {
    return mVariantStatistics.getStatistics();
  }

  /**
   * Prepare an string with an output line for one mutation
   * @param h1 writer for the first haploid set
   * @param h2 writer for the second haploid set (may be null)
   * @param seqName sequence name
   * @param m mutation
   * @param prevRefNt the reference nucleotide before the mutation
   * @param verbose give extra debug information
   * @param newAndOld strings with old bases or insert or deletions and new mutated bases or ..
   * @return string for mutation
   * @throws IOException when write fails
   */
  VcfRecord mutationLine(SdfWriter h1, SdfWriter h2, String seqName, Mutation m, char prevRefNt, boolean verbose, String... newAndOld) throws IOException {
    if ((m == null)
        || (seqName == null)
        || (newAndOld == null) || (newAndOld.length > 3) || (newAndOld.length < 2)
        || (m.mHeterozygous && (newAndOld.length < 2 || newAndOld.length > 3))
        || (!m.mHeterozygous && newAndOld.length != 2)) {
      throw new IllegalArgumentException("Unpossible");
    }
    for (final String str : newAndOld) {
      assert str.length() > 0;
    }

    if (newAndOld[1].length() == 0) {
      throw new IllegalArgumentException("Programmerror; first mutation must be given: " + StringUtils.LS + "Mutation: " + m.makeString());
    }
    if (m.mHeterozygous && newAndOld[1].equals(newAndOld[2])) {
      throw new IllegalArgumentException("Programmerror; wrong heterozygous mutation: " + m.makeString());
    }

    final String name = m.mPloidy == Ploidy.HAPLOID ? newAndOld[1] : newAndOld[1] + ":" + (m.mHeterozygous ? newAndOld[2] : newAndOld[1]);
    final VariantSample sample = new VariantSample(m.mPloidy, name.replace("i", ""), false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);

    writeToOutput(h1, newAndOld[1]);
    if (newAndOld.length > 2) {
      writeToOutput(h2, newAndOld[2]);
    } else if (!m.mHeterozygous) {
      writeToOutput(h2, newAndOld[1]);
    }
    final String refAllele = newAndOld[0].replace("i", "");
    final VariantLocus locus = new VariantLocus(seqName, m.mPos, m.mPos + refAllele.length(), refAllele, prevRefNt);
    final Variant v = new Variant(locus, sample);
    final VcfRecord rec = mFormatter.makeVcfRecord(v);
    if (mHeader != null) {
      mVariantStatistics.tallyVariant(mHeader, rec);
    }
    return rec;
  }

  private void writeToOutput(SdfWriter output, String mut) throws IOException {
    if (output != null) {
      final byte[] buf = new byte[mut.length()];
      int buflength = 0;
      for (int i = 0; i < mut.length(); i++) {
        final char dna = mut.charAt(i);
        if (dna != 'i' && dna != 'D') {
          buf[buflength] = (byte) DNA.valueOf(mut.charAt(i)).ordinal();
          buflength++;
        }
      }
      if (buflength > 0) {
        output.write(buf, null, buflength);
      }
    }
  }

  /** @return a detailed human readable representation of this object. It is intended that this show internal details of the object structure that may be relevant to an implementor/debugger but not to a user. */
  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("");
    if (mPriors != null) {
      sb.append(mPriors.toString()).append(StringUtils.LS);
      sb.append("====").append(StringUtils.LS);
      sb.append(mMutationGenerator.toString()).append(StringUtils.LS);
    }
    return sb.toString();
  }
}
