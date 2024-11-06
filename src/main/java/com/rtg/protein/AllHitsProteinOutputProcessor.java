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
        final ArrayList<File> alignmentFileList = new ArrayList<>();
        for (final ProteinOutputProcessor child : mChildren) {
          alignmentFileList.add(child.mOutFile);
        }
        final File[] files = alignmentFileList.toArray(new File[0]);
        if (quickcat) {
          FileUtils.copyRaw(mOutFile, files);
        } else {
          try (OutputStream out = FileUtils.createOutputStream(mOutFile)) {
            TsvUtils.tsvCat(gztemp, out, files);
          }
        }
        cleanup(mParams, alignmentFileList);
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
