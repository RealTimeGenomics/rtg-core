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
package com.rtg.metagenomics;

import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.usage.UsageMetric;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.SamRecordPopulator;

import net.sf.samtools.SAMRecord;

/**
 * Count mappings to a target genome versus all other genomes.
 */
@TestClass("com.rtg.metagenomics.SpeciesTargetCliTest")
class SpeciesTargetTask extends SpeciesTask {
  
  private final int mTargetSpeciesId;

  public SpeciesTargetTask(final SpeciesParams params, final OutputStream out, final UsageMetric usageMetric, final int targetId) {
    super(params, out, usageMetric);
    mTargetSpeciesId = targetId;
  }
  
  private final HashMap<String, long[]> mHits = new HashMap<>();

  @Override
  protected void accumulateMappings(final SequencesReader sr) throws IOException {
    Diagnostic.progress("Starting to read SAM records.");
    long usageStats = 0L;
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(mParams.mapped(), mParams.ioThreads(), new SingletonPopulatorFactory<>(new SamRecordPopulator()), mParams.filterParams(), SamUtils.getUberHeader(mParams.mapped()))) {
      while (it.hasNext()) {
        usageStats++;
        final SAMRecord rec = it.next();
        final String readId = rec.getReadName();
        final String sequenceName = rec.getReferenceName();
        final Integer taxonId = mSequenceMap.get(sequenceName);
        if (taxonId == null) {
          // something wrong - maybe mappings against different reference ?
          Diagnostic.warning("Could not find taxon ID for sequence: " + sequenceName);
        } else {
          long[] hitList = mHits.get(readId);
          if (hitList == null) {
            hitList = new long[2];
            mHits.put(readId, hitList);
          }
          final Integer speciesId = mSpeciesMap.id(taxonId);
          hitList[speciesId == mTargetSpeciesId ? 0 : 1]++;
        }
      }
    }
    mUsageMetric.setMetric(usageStats);
    Diagnostic.progress("Finished reading SAM records.");
  }

}
