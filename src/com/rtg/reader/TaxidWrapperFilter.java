/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import java.io.IOException;
import java.util.Collection;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.taxonomy.TaxonomyUtils;
import com.rtg.util.MultiMap;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Facilitates transferring sequences by (short) name from a SdfReaderWrapper to a WriterWrapper.
 */
@TestClass("com.rtg.reader.Sdf2FastaTest")
class TaxidWrapperFilter extends WrapperFilter {

  private final MultiMap<Integer, Long> mTaxToSeqId;

  TaxidWrapperFilter(SdfReaderWrapper reader, WriterWrapper writer) throws IOException {
    super(reader, writer);
    if (reader.isPaired()) {
      throw new NoTalkbackSlimException("Taxonomy is not supported for paired-end SDFs");
    }
    if (!TaxonomyUtils.hasTaxonomyInfo(reader.single())) {
      throw new NoTalkbackSlimException("The supplied SDF does not contain taxonomy information");
    }
    mTaxToSeqId = TaxonomyUtils.loadTaxonomyIdMapping(reader.single());
  }

  @Override
  protected void warnInvalidSequence(String seqid) {
    if (mWarnCount < 5) {
      Diagnostic.warning("No sequence data for taxonomy id " + seqid);
      mWarnCount++;
      if (mWarnCount == 5) {
        Diagnostic.warning("(Only the first 5 messages shown.)");
      }
    } else {
      Diagnostic.userLog("No sequence data for taxonomy id " + seqid);
    }
  }

  /**
   * Transfer an interpreted sequence or set of sequences from the reader to the writer.
   * This implementation interprets the specifier as a short sequence name.
   * @param seqRange the sequence name
   * @throws IOException if there was a problem during writing
   */
  @Override
  protected void transfer(String seqRange) throws IOException {
    final Collection<Long> ids = mTaxToSeqId.get(Integer.parseInt(seqRange));
    if (ids != null) {
      for (Long id : ids) {
        transfer(id);
      }
    } else {
      warnInvalidSequence(seqRange);
    }
  }
}
