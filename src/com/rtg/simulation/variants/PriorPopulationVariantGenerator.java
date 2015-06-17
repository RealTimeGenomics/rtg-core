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
package com.rtg.simulation.variants;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.NavigableSet;

import com.rtg.reader.SequencesReader;
import com.rtg.simulation.variants.PopulationMutatorPriors.VariantType;
import com.rtg.util.PortableRandom;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.array.longindex.LongChunks;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.SequenceIdLocus;
import com.rtg.util.intervals.SequenceIdLocusSimple;

/**
 * Generate population variants with types and distribution according to genome priors
 */
public class PriorPopulationVariantGenerator extends PopulationVariantGenerator {

  private static final class RandomPositionGenerator implements VariantPositionGenerator {

    private final SequencesReader mReader;
    private final PortableRandom mRandom;
    private final int[] mSequenceLengths;
    private final long mTotalLength;

    private RandomPositionGenerator(SequencesReader reader, PortableRandom random) throws IOException {
      mReader = reader;
      mRandom = random;
      mSequenceLengths = mReader.sequenceLengths(0, (int) mReader.numberSequences());
      mTotalLength = mReader.totalLength();
    }

    @Override
    public SequenceIdLocus nextVariantPosition() {
      long pos = (long) (mRandom.nextDouble() * mTotalLength);
      for (int i = 0; i < mSequenceLengths.length; i++) {
        if (pos < mSequenceLengths[i]) {
          return new SequenceIdLocusSimple(i, (int) pos);
        }
        pos -= mSequenceLengths[i];
      }
      return null;
    }
  }

  interface AlleleFrequencyChooser {
    double chooseAltFrequency();
  }

  static final class FixedAlleleFrequencyChooser implements AlleleFrequencyChooser {
    private final double mFreq;
    public FixedAlleleFrequencyChooser(double freq) {
      mFreq = freq;
    }
    @Override
    public double chooseAltFrequency() {
      return mFreq;
    }
  }

  static final class TableAlleleFrequencyChooser implements AlleleFrequencyChooser {

    private final PortableRandom mRandom;
    private final long[] mCumCounts;
    private final double[] mAltFreqs;

    // bias: tendency to emphasise alt alleles
    // from -1 (make alt alleles rarer) to 1 (make alt alleles much more common). 0 is no bias (effectively unaltered table).
    static AlleleFrequencyChooser make(PortableRandom random, double bias) {
      if (bias < -1 || bias > 1) {
        throw new IllegalArgumentException("Bias must be between -1 and 1");
      }
      try {
        try (BufferedReader r = new BufferedReader(new InputStreamReader(Resources.getResourceAsStream("com/rtg/simulation/variants/alt.allele.frequency.counts.txt")))) {
          String line;
          final LongChunks cumcounts = new LongChunks(0);
          final IntChunks altfreqs = new IntChunks(0);
          final double biasFactor = (bias + 1) / 2;
          double lastfreq = 0;
          long cumcount = 0;
          while ((line = r.readLine()) != null) {
            line = line.trim();
            if ((line.length() == 0) || (line.charAt(0) == '#')) {
              continue;
            }
            final String[] parts = StringUtils.split(line, ' ');
            if (parts.length != 2) {
              throw new RuntimeException("Malformed allele frequency priors");
            }

            final double freq = Double.valueOf(parts[0]);
            assert freq >= lastfreq : "Out-of-order allele frequency in priors " + freq;
            assert freq >= 0 : "Illegal allele frequency in priors " + freq;
            assert freq <= 1 : "Illegal allele frequency in priors " + freq;

            double count = Double.valueOf(parts[1]);
            count = biasFactor * count * freq + (1.0 - biasFactor) * count * (1 - freq);
            cumcount += count;
            assert count >= 0 : "Illegal negative allele frequency count in priors " + count;
            assert cumcount >= 0 : "Overflow in allele frequency prior counts";

            altfreqs.append((int) (1000 * freq)); // 3dp of precision in the freq
            cumcounts.append(cumcount);
            lastfreq = freq;
          }
          final long[] cumCountsArr = new long[(int) cumcounts.length()];
          final double[] altFreqsArr = new double[(int) cumcounts.length()];
          for (int i = 0; i < cumcounts.length(); i++) {
            cumCountsArr[i] = cumcounts.get(i);
            altFreqsArr[i] = (double) altfreqs.get(i) / 1000;
          }
          return new TableAlleleFrequencyChooser(random, cumCountsArr, altFreqsArr);
        }
      } catch (IOException e) {
        throw new RuntimeException("Could not load allele frequency priors");
      }
    }

    TableAlleleFrequencyChooser(PortableRandom random, long[] cumCounts, double[] altFreqs) {
      mRandom = random;
      mCumCounts = cumCounts;
      mAltFreqs = altFreqs;
    }

    @Override
    public double chooseAltFrequency() {
      final long pos = (long) (mRandom.nextDouble() * mCumCounts[mCumCounts.length - 1]);
      int idx = Arrays.binarySearch(mCumCounts, pos);
      if (idx < 0) {
        idx = -(idx + 1);
      }
      return mAltFreqs[idx];
    }
  }

  private final SequencesReader mReader;
  private final PortableRandom mRandom;
  private final PopulationMutatorPriors mPriors;
  private final AlleleFrequencyChooser mFreqChooser;
  private final byte[] mTemplate;
  private final long mTargetVariantCount;
  private double mAcceptedCount = 0;


  /**
   * @param reader reference sequences
   * @param priors population mutator priors
   * @param random source of randomness
   * @param bias bias frequencies toward alt alleles. 0 is no bias. 1 is fully biased
   * @throws IOException if an I/O error occurs.
   */
  public PriorPopulationVariantGenerator(SequencesReader reader, PopulationMutatorPriors priors, PortableRandom random, double bias) throws IOException {
    this(reader, priors, random, bias, (int) (reader.totalLength() * priors.rate()));
  }

  /**
   * @param reader reference sequences
   * @param priors population mutator priors
   * @param random source of randomness
   * @param bias bias frequencies toward alt alleles. 0 is no bias. 1 is fully biased
   * @param targetVariants the number of variants a typical member of the population should end up with
   * @throws IOException if an I/O error occurs.
   */
  public PriorPopulationVariantGenerator(SequencesReader reader, PopulationMutatorPriors priors, PortableRandom random, double bias, int targetVariants) throws IOException {
    this(reader, priors, random, TableAlleleFrequencyChooser.make(random, bias), targetVariants);
  }

  /**
   * @param reader reference sequences
   * @param priors population mutator priors
   * @param random source of randomness
   * @param freqChooser chooses allele frequencies
   * @param targetVariants the number of variants a typical member of the population should end up with
   * @throws IOException if an I/O error occurs.
   */
  public PriorPopulationVariantGenerator(SequencesReader reader, PopulationMutatorPriors priors, PortableRandom random, AlleleFrequencyChooser freqChooser, int targetVariants) throws IOException {
    super(new RandomPositionGenerator(reader, random));
    mReader = reader;
    mPriors = priors;
    mRandom = random;
    mFreqChooser = freqChooser;
    mTargetVariantCount = targetVariants;
    mTemplate = new byte[(int) reader.maxLength()];
    Diagnostic.developerLog("Generating " + mTargetVariantCount + " expected variants per genome");
  }

  @Override
  protected boolean checkValid(PopulationVariant var, NavigableSet<PopulationVariant> set) {
    final boolean accepted = super.checkValid(var, set);
    if (accepted) {
      double sum = 0;
      for (int i = 0; i < var.mDistribution.length; i++) {
        sum += var.mDistribution[i];
      }
      mAcceptedCount += sum; // 1 - p(ref)
    }
    return accepted;
  }

  @Override
  protected boolean needMoreVariants() {
    // Once sum of AF of accepted variants >= mTargetVariantCount, then stop
    return mAcceptedCount < mTargetVariantCount;
  }

  @Override
  PopulationVariant nextPopulationVariant() throws IOException {
    final SequenceIdLocus vp = mVariantPositionGenerator.nextVariantPosition();
    if (vp == null) {
      return null;
    }
    //System.err.println(vp.getReferenceId() + ":" + position);

    // Now choose the variant according to the priors
    final PopulationVariant popVar = new PopulationVariant(vp);

    // Currently only simulate biallelic variants
    popVar.mAlleles = new byte[1][];
    popVar.mDistribution = new double[] {mFreqChooser.chooseAltFrequency()};

    fillAlleles(popVar);
    return popVar;
  }

  private void fillAlleles(final PopulationVariant popVar) throws IOException {
    VariantType type = mPriors.chooseType(mRandom);
    int variantLength = mPriors.chooseLength(mRandom, type);
    final int position = popVar.getStart();
    if (position + variantLength >= mReader.length(popVar.getSequenceId())) {
      type = VariantType.SNP;
      variantLength = 1;
    }
    switch (type) {
      case INSDEL: // XXX Temporarily treat insdel as SNP. Need to implement this
        //System.err.println("Insdel wants length " + variantLength);
      case SNP:
        mReader.read(popVar.getSequenceId(), mTemplate, position, 1);
        popVar.mRef = Arrays.copyOf(mTemplate, 1);
        popVar.mAlleles[0] = mPriors.chooseAltSnp(mRandom, popVar.mRef[0]);
        break;
      case MNP:
        assert variantLength > 1 : "MNP needs variant length greater than 1";
        mReader.read(popVar.getSequenceId(), mTemplate, position, variantLength);
        popVar.mRef = Arrays.copyOf(mTemplate, variantLength);
        popVar.mAlleles[0] = mPriors.chooseAltMnp(mRandom, popVar.mRef);
        //System.err.println("MNP " + DnaUtils.bytesToSequenceIncCG(popVar.mRef) + " -> " + DnaUtils.bytesToSequenceIncCG(popVar.mAlleles[0]));
        break;
      case DELETE:
        assert variantLength > 0 : "Deletions needs variant length greater than 1";
        mReader.read(popVar.getSequenceId(), mTemplate, position, variantLength + 1); // Need an anchor base
        popVar.mRef = Arrays.copyOf(mTemplate, variantLength + 1);
        popVar.mAlleles[0] = new byte[] {popVar.mRef[0]};
        //System.err.println("Delete " + DnaUtils.bytesToSequenceIncCG(popVar.mRef) + " -> " + DnaUtils.bytesToSequenceIncCG(popVar.mAlleles[0]));
        break;
      case INSERT:
        assert variantLength > 0 : "Insertions needs variant length greater than 1";
        mReader.read(popVar.getSequenceId(), mTemplate, position, variantLength + 1); // Need anchor base
        popVar.mRef = Arrays.copyOf(mTemplate, 1);

        popVar.mAlleles[0] = mPriors.chooseAltMnp(mRandom, Arrays.copyOf(mTemplate, variantLength + 1)); //XXX Improve this
        popVar.mAlleles[0][0] = popVar.mRef[0];
        //System.err.println("Insert " + DnaUtils.bytesToSequenceIncCG(popVar.mRef) + " -> " + DnaUtils.bytesToSequenceIncCG(popVar.mAlleles[0]));
        break;
/*
      case INSDEL:
        mReader.read(vp.getReferenceId(), mTemplate, position, variantLength);
        break;
*/
      default:
        throw new IllegalArgumentException("Invalid mutation distribution");
    }
  }
}
