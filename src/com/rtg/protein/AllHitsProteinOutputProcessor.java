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
package com.rtg.protein;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.TreeSet;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.MapStatistics;
import com.rtg.ngs.NgsParams;
import com.rtg.util.TsvUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.GzipAsynchOutputStream;

/**
 * Thread safe implementation of All Hits output processor for protein
 * world, taking care of duplicate hit removal on a per-template
 * basis.
 *
 */
@TestClass(value = {"com.rtg.protein.ProteinOutputProcessorTest"})
public class AllHitsProteinOutputProcessor extends ProteinOutputProcessor {

  private final TreeSet<ProteinAlignmentResult> mTemplateHits = new TreeSet<>();
  private int mBufferDistance = 1;
  private ProteinAlignmentResult mLastWrittenRecord = null;

  /**
   * Construct a new {@link AllHitsProteinOutputProcessor}
   * @param params {@link NgsParams} object
   * @param statistics collector of statistics.
   * @throws IOException if error
   */
  public AllHitsProteinOutputProcessor(NgsParams params, MapStatistics statistics) throws IOException {
    super(params, statistics);
    mBufferDistance = (int) mParams.buildFirstParams().reader().maxLength(); // Will be plenty, as mBufferDistance units are amino acids, but build is nucleotides
    assert params.buildFirstParams().numberSequences() <= Integer.MAX_VALUE;
  }

  /**
   * Create a new {@link AllHitsProteinOutputProcessor} for use as a child processor during multi-core processing.
   */
  private AllHitsProteinOutputProcessor(NgsParams params, AllHitsProteinOutputProcessor master, int childId, SharedStatusCollector collector, MapStatistics statistics) throws IOException {
    super(params, master, childId, false, collector, statistics);
    mBufferDistance = (int) mParams.buildFirstParams().reader().maxLength(); // Will be plenty, as mBufferDistance units are amino acids, but build is nucleotides
  }

  @Override
  public void finish() throws IOException {
    if ((this == mMaster) && mChildren.size() > 0) {
      closeChildren();
      if (mParams.outputParams().mergeAlignmentResults()) {
        close();  // needs to release open file, because silly master opens output file even when multithreaded (due to subclasses)
        final boolean gztemp = mParams.outputParams().isCompressOutput();
        final boolean quickcat = mParams.outputParams().isCompressOutput() == gztemp; // If both, we can avoid decompress/recompress during final cat
        OutputStream out = null;
        final ArrayList<File> alignmentFileList = new ArrayList<>();
        try {
          for (final ProteinOutputProcessor child : mChildren) {
            alignmentFileList.add(child.mOutFile);
          }
          final File[] files = alignmentFileList.toArray(new File[alignmentFileList.size()]);
          if (quickcat) {
            FileUtils.catInSync(mOutFile, false, files);
          } else {
            if (mParams.outputParams().isCompressOutput()) {
              out = new GzipAsynchOutputStream(mOutFile);
            } else {
              out = FileUtils.createOutputStream(mOutFile, false, false);
            }
            TsvUtils.tsvCat(gztemp, out, files);
          }
        } finally {
          if (out != null) {
            out.flush();
            out.close();
          }
          cleanup(mParams, alignmentFileList);
        }
      }
    } else {
      flushTemplateHits();
      close();
    }
    if (this == mMaster) {
      writeUnmapped();
      mSharedStatusCollector.calculateStatistics();
    }
  }

  @Override
  protected void nextTemplateId(long templateId) throws IOException {
    flushTemplateHits();
    super.nextTemplateId(templateId);
  }

  /** Write out the hits for the current template, after removing duplicates */
  private void flushTemplateHits() throws IOException {
    //System.out.println("Flushing " + mTemplateHits.size() + " template hits for " + mCurrentTemplateId);
    for (ProteinAlignmentResult hit : mTemplateHits) {
      super.writeResult(hit);
    }
    mTemplateHits.clear();
    mLastWrittenRecord = null;
  }

  /**
   * Rather than writing the result directly, this implementation
   * buffers all results for the current template, to be written later
   * after duplicate removal.
   */
  @Override
  void writeResult(final ProteinAlignmentResult res) throws IOException {
    if (mLastWrittenRecord == null
        || mLastWrittenRecord.getAlignmentStart() < res.getAlignmentStart()) { // Must be < not <= here to ensure we can still detect duplicates.
      mTemplateHits.add(res);
      final int lastPos = mTemplateHits.last().getAlignmentStart();
      while (lastPos - mTemplateHits.first().getAlignmentStart() > mBufferDistance) {
        final ProteinAlignmentResult firstRecord = mTemplateHits.first();
        mTemplateHits.remove(firstRecord);
        super.writeResult(firstRecord);
        mLastWrittenRecord = firstRecord;
      }
    } else {
      throw new IllegalStateException("all-hits writer buffer distance (" + mBufferDistance + " bases, " + mTemplateHits.size() + " records buffered)"
                                      + " too small by " + (mLastWrittenRecord.getAlignmentStart() - res.getAlignmentStart()));
    }
    //super.writeResult(res);
  }

  private void cleanup(final NgsParams params, final ArrayList<File> alignmentFileList) {
    if (!params.outputParams().keepIntermediate()) {
      for (final File f : alignmentFileList) {
        if (f.delete()) {
          Diagnostic.userLog("Deleted file: " + f.getPath());
        } else {
          Diagnostic.userLog("Cannot delete file: " + f.getPath());
        }
      }
    }
  }

  @Override
  public OutputProcessor threadClone(HashingRegion region) throws IOException {
    if (region != HashingRegion.NONE) {
      throw new UnsupportedOperationException();
    }
    final int current = mChildren.size();
    final AllHitsProteinOutputProcessor p = new AllHitsProteinOutputProcessor(super.getParams(), this, current, mSharedStatusCollector, mStatistics);
    mChildren.add(p);
    return p;
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      flushTemplateHits();
    } finally {
      close();
    }
  }
}
