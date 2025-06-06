/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.ngs;

import java.util.Collection;
import java.util.Collections;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.alignment.AlignerMode;
import com.rtg.alignment.EditDistanceFactory;
import com.rtg.index.FixedRepeatFrequencyFilterMethod;
import com.rtg.index.IndexFilterMethod;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.ModuleParams.ModuleParamsBuilder;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.protein.ProteinReadIndexer;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.test.params.BuilderNotSet;

/**
 * A builder class for <code>NgsParams</code>.
 */
@TestClass(value = "com.rtg.ngs.NgsParamsTest")
public class NgsParamsBuilder extends ModuleParamsBuilder<NgsParamsBuilder> {

  int mNumberThreads = 1;
  int mThreadMultiplier = HashingRegion.DEFAULT_THREAD_MULTIPLIER;
  Integer mMinFragmentLength = null;
  Integer mMaxFragmentLength = null;
  //  Integer mExpectedInsertSize = null;
  int mReadFreqThreshold = 65535;
  int mStepSize = -1;
  boolean mUseLongReadMapping = false;
  ProteinScoringMatrix mProteinScoringMatrix = null;
  boolean mEnableProteinReadCache = false;
  ISequenceParams mBuildFirstParams = null;
  ISequenceParams mBuildSecondParams = null;
  ISequenceParams mSearchParams = null;
  NgsOutputParams mOutputParams = null;
  Collection<ListenerType> mListeners = Collections.singleton(ListenerType.NULL);
  boolean mCompressHashes = true;
  int mIntSetWindow = 1;
  Integer mMinHits = null;
  boolean mLegacyCigars = false;
  boolean mUseTopRandom = false;
  long mMapXMinReadLength = -1;
  NgsMaskParams mMaskParams = new NgsMaskParamsGeneral(MapFlags.DEFAULT_WORD_SIZE, 1, 0, 1);
  int mMapXMetaChunkSize = ProteinReadIndexer.DEFAULT_META_CHUNK_LENGTH;
  int mMapXMetaChunkOverlap = ProteinReadIndexer.DEFAULT_META_CHUNK_OVERLAP;
  MachineOrientation mPairOrientation = MachineOrientation.ANY;
  boolean mParallelUnmatedProcessing = false;
  int mGapOpenPenalty = EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY;
  int mGapExtendPenalty = EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY;
  int mSubstitutionPenalty = EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY;
  int mUnknownsPenalty = EditDistanceFactory.DEFAULT_UNKNOWNS_PENALTY;
  int mIndelSoftClipDistance = 0;
  int mMismatchSoftClipDistance = 0;
  int mMinMatches = 0;
  MaxShiftFactor mAlignerBandWidthFactor = new MaxShiftFactor(0.5);
  AlignerMode mAlignerMode = AlignerMode.AUTO;
  String mSingleIndelPenalties = EditDistanceFactory.DEFAULT_SINGLE_INDEL_TABLE;
  IndexFilterMethod mIndexFilter = new FixedRepeatFrequencyFilterMethod(1000);


  @Override
  protected NgsParamsBuilder self() {
    return this;
  }

  /**
   * Sets whether long read mapping should be used.
   *
   * @param useLong true if long read mapping should be used.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder useLongReadMapping(final boolean useLong) {
    mUseLongReadMapping = useLong;
    return self();
  }

  /**
   * Sets the number of threads to use.
   *
   * @param numThreads the number of threads to employ. The default is 1.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder numberThreads(final int numThreads) {
    mNumberThreads = numThreads;
    return self();
  }

  /**
   * Sets the multiplier factor when dividing work into chunks among multiple threads.
   *
   * @param multiplier the number of chunks of work per thread.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder threadMultiplier(final int multiplier) {
    mThreadMultiplier = multiplier;
    return self();
  }

  /**
   * Sets the minimum fragment length for paired reads to mate.
   *
   * @param size the minimum fragment length. The default is 0.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder minFragmentLength(final int size) {
    mMinFragmentLength = size;
    return self();
  }

  /**
   * Sets the maximum fragment length for paired reads to mate.
   *
   * @param size the maximum fragment length. The default is 0.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder maxFragmentLength(final int size) {
    mMaxFragmentLength = size;
    return self();
  }

  /**
   * Sets the maximum number of hits a read may have before it is discarded.
   *
   * @param threshold the threshold. The default is 65535.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder readFreqThreshold(final int threshold) {
    mReadFreqThreshold = threshold;
    return self();
  }

  /**
   * Sets the step size.
   *
   * @param step the step size. The default is 1.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder stepSize(final int step) {
    mStepSize = step;
    return self();
  }

  /**
   * Sets the build parameters.
   *
   * @param params the parameters used during build.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder buildFirstParams(final ISequenceParams params) {
    mBuildFirstParams = params;
    return self();
  }

  /**
   * @return the build and create parameters.
   */
  @BuilderNotSet
  public ISequenceParams buildFirstParams() {
    return mBuildFirstParams;
  }

  /**
   * Sets the build parameters used for the right side when in paired-end mode.
   *
   * @param params the parameters used during build.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder buildSecondParams(final ISequenceParams params) {
    mBuildSecondParams = params;
    return self();
  }

  /**
   * @return the build and create parameters.
   */
  @BuilderNotSet
  public ISequenceParams buildSecondParams() {
    return mBuildSecondParams;
  }

  /**
   * Sets the search parameters.
   *
   * @param params the parameters used during search.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder searchParams(final ISequenceParams params) {
    mSearchParams = params;
    return self();
  }

  /**
   * @return the search parameters.
   */
  @BuilderNotSet
  public ISequenceParams searchParams() {
    return mSearchParams;
  }

  /**
   * Sets the output parameters.
   *
   * @param params the parameters used for output.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder outputParams(final NgsOutputParams params) {
    mOutputParams = params;
    return self();
  }

  /**
   * @param mask parameters for the mask used in index/searching
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder maskParams(NgsMaskParams mask) {
    mMaskParams = mask;
    return self();
  }
  /**
   * @param windowSelector select which split window class to use.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder maskParams(final String windowSelector) {
    mMaskParams = new NgsMaskParamsExplicit(windowSelector);
    return self();
  }

  /**
   * Sets listeners.
   *
   * @param listeners a collection of listeners.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder listeners(final Collection<ListenerType> listeners) {
    mListeners = listeners;
    return self();
  }

  /**
   * Sets protein scoring matrix.
   *
   * @param matrix the protein scoring matrix.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder proteinScoringMatrix(final ProteinScoringMatrix matrix) {
    mProteinScoringMatrix = matrix;
    return self();
  }

  /**
   * Sets whether protein read cache should be used.
   *
   * @param value true if protein read cache should be used.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder enableProteinReadCache(boolean value) {
    mEnableProteinReadCache = value;
    return self();
  }

  /**
   * Sets compress hashes flag
   * @param compressHashes the compress hashes flag
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder compressHashes(final boolean compressHashes) {
    mCompressHashes = compressHashes;
    return self();
  }

  /**
   * Sets <code>IntSet</code> window for hashes
   * @param window window within which hashes does not call call method
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder intsetWindow(final int window) {
    mIntSetWindow = window;
    return self();
  }

  /**
   * If not using gapped output for long reads this parameter
   * determines how many hits to a region a read needs before it
   * is allowed to read output processing.
   * @param minHits the number of hits
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder minHits(Integer minHits) {
    mMinHits = minHits;
    return self();
  }

  /**
   * Whether to output legacy cigar format
   * @param legacy whether to output legacy cigar format
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder legacyCigars(boolean legacy) {
    mLegacyCigars = legacy;
    return self();
  }

  /**
   * @param useTopRandom if true use top random output.
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder useTopRandom(final boolean useTopRandom) {
    mUseTopRandom = useTopRandom;
    return self();
  }

  /**
   * @param minLength if true treat N's as mismatches in alignments
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder mapXMinLength(final long minLength) {
    mMapXMinReadLength = minLength;
    return self();
  }

  /**
   * @param size the size at which meta chunks are created
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder mapXMetaChunkSize(int size) {
    mMapXMetaChunkSize = size;
    return self();
  }

  /**
   * @param overlap how much overlap to have between adjacent meta chunks
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder mapXMetaChunkOverlap(int overlap) {
    mMapXMetaChunkOverlap = overlap;
    return self();
  }

  /**
   * @param o correct orientation for mating
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder pairOrientation(MachineOrientation o) {
    mPairOrientation = o;
    return self();
  }

  /**
   * @param value set true if parallelization of the mated+unmated processing stage is desired
   * @return this builder, so calls can be chained.
   */
  public NgsParamsBuilder parallelUnmatedProcessing(boolean value) {
    mParallelUnmatedProcessing = value;
    return self();
  }

  /**
   * @param value the penalty for a gap open during alignment
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder gapOpenPenalty(int value) {
    mGapOpenPenalty = value;
    return self();
  }

  /**
   * @param value the penalty for a gap extension during alignment
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder gapExtendPenalty(int value) {
    mGapExtendPenalty = value;
    return self();
  }

  /**
   * @param value the penalty for a substitution during alignment
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder substitutionPenalty(int value) {
    mSubstitutionPenalty = value;
    return self();
  }

  /**
   * @param value the penalty for an unknown nucleotide during alignment
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder unknownsPenalty(int value) {
    mUnknownsPenalty = value;
    return self();
  }

  /**
   * @param value soft clip alignments if an indel occurs VALUE bp from either end
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder indelSoftClipDistance(int value) {
    mIndelSoftClipDistance = value;
    return self();
  }

  /**
   * @param value soft clip alignments if a mismatch occurs VALUE bp from either end
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder mismatchSoftClipDistance(int value) {
    mMismatchSoftClipDistance = value;
    return self();
  }

  /**
   * @param value the minimum number of matches that must be present after soft clipping for an alignment to be retained.
   * @return this builder, so calls can be chained
   */
  public NgsParamsBuilder minMatches(int value) {
    mMinMatches = value;
    return self();
  }

  /**
   * @param value the aligner max shift band factor
   * @return this build, so calls can be chained
   */
  public NgsParamsBuilder alignerBandWidthFactor(MaxShiftFactor value) {
    mAlignerBandWidthFactor = value;
    return self();
  }
  /**
   * @param value the aligner chain to use
   * @return this build, so calls can be chained
   */
  public NgsParamsBuilder alignerMode(AlignerMode value) {
    mAlignerMode = value;
    return self();
  }
  /**
   * @param value the name of the properties resource to load single indel table from
   * @return this build, so calls can be chained
   */
  public NgsParamsBuilder singleIndelPenalties(String value) {
    mSingleIndelPenalties = value;
    return self();
  }

  /**
   * @param filter handle that is used to filter (normally repetitive) hashes from index
   * @return this build, so calls can be chained
   */
  public NgsParamsBuilder indexFilter(IndexFilterMethod filter) {
    mIndexFilter = filter;
    return self();
  }

  /**
   * Creates a NgsParams using the current builder
   * configuration.
   * @return the new NgsParams
   */
  public NgsParams create() {
    return new NgsParams(this);
  }
}
