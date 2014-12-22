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
package com.rtg.ngs.tempstage;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.alignment.EditDistance;
import com.rtg.alignment.EditDistanceFactory;
import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Base class used to derive temp file writers.
 */
@TestClass(value = {"com.rtg.ngs.tempstage.PairedTempFileWriterImplTest", "com.rtg.ngs.tempstage.SingleEndTempFileWriterTest"})
public abstract class AbstractTempFileWriter implements Closeable {

  /**
   * Turning this on allows for tighter max score bounds as we find good hits for a read, which means aligners can terminate earlier.
   * Unfortunately Hop Step (and probably Seeded) are a bit too heuristic and cause non-determinism issues when they find 'worse' hits
   * at a higher threshold, but when given a lower threshold, pass on to another aligner which can find a better alignment.
   * In other words, if we can make a more accurate aligner chain, turning this on would be beneficial.
   */
  private static final boolean USE_BLOCKERS_FOR_EARLY_TERM = false; //Boolean.valueOf(System.getProperty("rtg.blockers-for-early-term", "false"));

  private static final boolean DUMP_ALIGNMENT_STATS = GlobalFlags.isSet(GlobalFlags.TEMP_FILES_DUMP_ALIGN_STATS_FLAG);
  protected final SharedResources mSharedResources;
  protected final ReadStatusListener mListener;
  protected SequencesReader mTemplateReader;
  protected SequencesReader mFirstReader;
  protected SequencesReader mSecondReader;
  protected final IntegerOrPercentage mMatedMaxMismatches;
  protected final IntegerOrPercentage mUnmatedMaxMismatches;
  protected HashingRegion mClipRegion = HashingRegion.NONE;
  protected byte[] mTemplate = null;
  protected int mTemplateOffset = 0;
  protected long mTemplateId = -1;
  protected EditDistance mMrEd;
  protected final int mTemplatePadding;
  protected final boolean mLegacy;
  protected final int mSubstitutionPenalty;
//  protected final HashingRegion mReadRange;

  private int mLastReadLength = -1;
  private int mLastReadLengthMaxShiftValue = -1;

  private final MaxShiftFactor mMaxShiftFactor;

  /**
   * Construct a new writer.
   * @param listener a listener that will receive notifications of mated pairs
   * @param sharedResources resources shared between multiple instances of this class
   * @param params {@link NgsParams} for current run
   */
  public AbstractTempFileWriter(ReadStatusListener listener, SharedResources sharedResources, final NgsParams params) {
    mSharedResources = sharedResources;
    mFirstReader = mSharedResources.firstReaderCopy();
    mSecondReader = mSharedResources.secondReaderCopy();
    mTemplateReader = mSharedResources.templateReaderCopy();
    mMatedMaxMismatches = params.outputParams().filter().matedMaxMismatches();
    mUnmatedMaxMismatches = params.outputParams().filter().unmatedMaxMismatches();
//    mReadRange = params.buildFirstParams().readerRestriction();

    mMaxShiftFactor = params.alignerBandWidthFactor();
    mMrEd = EditDistanceFactory.createEditDistance(params, mFirstReader, mSecondReader);
    mListener = listener;
    mTemplatePadding = params.calculateThreadPadding();

    mLegacy = params.legacyCigars();
    mSubstitutionPenalty = params.substitutionPenalty();
  }

  /**
   * Set the region beyond which this thread should not report results
   * @param region the region to which results must be restricted.
   */
  public void setClipRegion(final HashingRegion region) {
    mClipRegion = region;
  }

  /**
   * Step to the specified template identifier.  The ids should be monotonically
   * increasing. <code>Long.MAX_VALUE</code> should be used to indicate no more sequences.
   *
   * @param templateId a <code>long</code> value
   * @throws IOException if an I/O Error occurs
   */
  public void nextTemplateId(final long templateId) throws IOException {
    if (templateId <= mTemplateId) {
      throw new IllegalArgumentException();
    }
    mTemplateId = templateId;
    if (templateId == Long.MAX_VALUE) {
      //      Diagnostic.developerLog("AbstractSamAlignmentWriter Long.MAX_VALUE");
      mTemplate = null;
    } else {

      mTemplateReader.seek(templateId);
      final int templateLength = mTemplateReader.currentLength();

      final int startpos = mClipRegion.getReferenceStart(templateId, mTemplatePadding);
      final int endpos = mClipRegion.getReferenceEnd(templateId, mTemplatePadding, templateLength);

      final int newLength = endpos - startpos;
      if ((newLength > 0) && (newLength < templateLength)) {
        Diagnostic.developerLog("Trimming template " + mTemplateId + " (0," + templateLength + ") to " + startpos + "," + endpos + " (length " + newLength + ") for region: " + mClipRegion + " using padding: " + mTemplatePadding);
        mTemplateOffset = startpos;
        //        Diagnostic.developerLog("AbstractSamAlignmentWriter Allocating: " + newLength + " bytes");
        mTemplate = new byte[newLength];
        if (mTemplateReader.readCurrent(mTemplate, startpos, newLength) != newLength) {
          throw new RuntimeException();
        }
      } else {
        mTemplateOffset = 0;
        //        Diagnostic.developerLog("AbstractSamAlignmentWriter Allocating (whole): " + templateLength + " bytes, startpos= " + startpos + " endpos= " + endpos + " newlength=" + newLength);
        mTemplate = new byte[templateLength];
        if (mTemplateReader.readCurrent(mTemplate) != templateLength) {
          throw new RuntimeException();
        }
      }
    }
  }

  protected int[] calculateEditDistance(byte[] read, int length, int start, boolean rc, IntegerOrPercentage maxMismatches, boolean left, int readId) {
    if (mLastReadLength != length) {
      mLastReadLength = length;
      mLastReadLengthMaxShiftValue = mMaxShiftFactor.calculateMaxShift(length);
    }

    final int score = maxMismatches.getValue(length) * mSubstitutionPenalty;
    final int leastScore = USE_BLOCKERS_FOR_EARLY_TERM ? Math.min(mSharedResources.getBlocker().getTerminationScore(readId), score) : score;

    return mMrEd.calculateEditDistance(read, length, mTemplate, start - mTemplateOffset, rc, leastScore, mLastReadLengthMaxShiftValue, left);
  }

  public MapQScoringReadBlocker getBlocker() {
    return mSharedResources.getBlocker();
  }

  protected SmartTempFileWriter createSmartWriter(final OutputStream out) {
    final int maxLength = mSecondReader == null ? (int) mFirstReader.maxLength() : Math.max((int) mFirstReader.maxLength(), (int) mSecondReader.maxLength());
    // Initial reported position can be out by maxshift for an individual record. Between records though the ordering could be out by read length + maxgap.

    final int maxMaxShift = mMaxShiftFactor.calculateMaxShift(maxLength);
//    final int maxShift = Math.max(mMaxShift, MaxShiftUtils.calculateDefaultMaxShift(maxLength)) * 1;

    return new SmartTempFileWriter(out, getRecordComparator(), maxLength * 2 + maxMaxShift);
  }

  abstract Comparator<BinaryTempFileRecord> getRecordComparator();


  /**
   * Closes resources used by this class.
   * @throws IOException if an IO exception occurs
   */
  @Override
  public void close() throws IOException {
    mSharedResources.close();

    if (mTemplateReader != null) {
      mTemplateReader.close();
      mTemplateReader = null;
    }
    if (mFirstReader != null) {
      mFirstReader.close();
      mFirstReader = null;
    }
    if (mSecondReader != null) {
      mSecondReader.close();
      mSecondReader = null;
    }
    if (DUMP_ALIGNMENT_STATS && mMrEd != null) {
      mMrEd.logStats();
    }
    mMrEd = null;

    Diagnostic.developerLog("AbstractSamAlignmentWriter.close() ");
    if (mTemplate != null) {
      Diagnostic.developerLog("Freeing: " + mTemplate.length + " bytes");
      mTemplate = null;
    }
  }

}
