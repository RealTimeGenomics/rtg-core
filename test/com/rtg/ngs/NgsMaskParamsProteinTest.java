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

import java.io.StringWriter;

import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallMock;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallMock;
import com.rtg.index.hash.ngs.protein.ProteinMask;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test for corresponding class
 */
public class NgsMaskParamsProteinTest extends TestCase {

  public void testMaskFactory() {
    Diagnostic.setLogStream();
    final NgsMaskParamsProtein params = new NgsMaskParamsProtein(9, 0, 0, 1);
    final ProteinMask mask = (ProteinMask) params.maskFactory(36).create(new ReadCallMock(new StringWriter()), new TemplateCallMock(new StringWriter()));
    assertEquals(11, mask.readLength());
  }

}
