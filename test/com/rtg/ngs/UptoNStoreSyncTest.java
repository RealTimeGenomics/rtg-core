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
package com.rtg.ngs;

import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class UptoNStoreSyncTest extends TopNImplementationTest {

  @Override
  protected UptoNStore getTopNImplementation(final int numReads, final int numTemplateSeqs, final int n, final int templateMaxLength) {
    return new UptoNStoreSync(new TopNImplementation(numReads, numTemplateSeqs, n, templateMaxLength, 50));
  }

  @Override
  public void testToString() {
    Diagnostic.setLogStream();
    final UptoNStore topn = getTopNImplementation(4, 17, 5, 1000);
    assertEquals("TopNImplementationSync", topn.toString());
  }

  public void testSyncId() {
    assertEquals(0, UptoNStoreSync.syncId(0));
    assertEquals(0, UptoNStoreSync.syncId(65536));
    assertEquals(1, UptoNStoreSync.syncId(1));
    assertEquals(1, UptoNStoreSync.syncId(65536 + 1));
  }
}
