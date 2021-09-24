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
package com.rtg.index.hash.ngs.instances;


import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallAccumulate;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallCheck;

import org.junit.Assert;

/**
 * Check tricky mask cases.
 */
public final class CheckMask {

  private CheckMask() { }

  static void check(final HashFunctionFactory factory, final String read, final String template) throws IOException {
    final ReadCallAccumulate rc = new ReadCallAccumulate();
    final TemplateCallCheck tcc = new TemplateCallCheck(rc.map());
    final NgsHashFunction hf = factory.create(rc, tcc);
    AbstractSplitTest.encode(hf, read);
    hf.readAll(0, false);
    hf.reset();
    AbstractSplitTest.encode(hf, template);
    hf.templateForward(0);
    hf.reset();
    tcc.done();
    Assert.assertTrue(tcc.mFoundDone);
  }
}
