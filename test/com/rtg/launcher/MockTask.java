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
package com.rtg.launcher;

import java.io.OutputStream;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;

/**
 */
public class MockTask extends ParamsTask<MockCliParams, NoStatistics> {

  MockTask(final MockCliParams params, final OutputStream defaultOutput) {
    super(params, defaultOutput, new NoStatistics(), null);
  }

  @Override
  protected void exec() {
    if (mParams.mRuntime) {
      throw new RuntimeException("RuntimeException thrown during execution of task");
    }
    if (mParams.mRuntimeSlim) {
      throw new SlimException(new RuntimeException("SlimException thrown during execution of task"));
    }
    if (mParams.mRuntimeErr) {
      Diagnostic.userLog("Write to err during task execution");
    }
    //new OutputStreamWriter(out).append("Mock did something" + LS);
  }
}
