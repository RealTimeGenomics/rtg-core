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
package com.rtg.util.cli;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.StringReader;

import com.reeltwo.spelling.Spelling;
import com.rtg.util.StringUtils;

import junit.framework.Assert;

/**
 * Separate class to keep Sharpen happy.
 *
 */
public final class CheckSpelling {

  static Spelling sSpelling = null;

  private CheckSpelling() { }

  static void setSpelling(final StringBuilder problems) throws IOException {
    sSpelling = new Spelling() {
        {
          includeStandardDictionaries();
          addCaseInsensitiveDictionary("com/rtg/util/cli/spell.insensitive");
          addCaseSensitiveDictionary("com/rtg/util/cli/spell.sensitive");
        }

        @Override
        protected void warning(final String source, final int lineNumber, final String msg) {
          problems.append("Flag: --").append(source).append(" ").append(msg).append(StringUtils.LS);
        }
      };
  }

  static void check(final String name, final String s) {
    try {
      try (LineNumberReader r = new LineNumberReader(new StringReader(s))) {
        sSpelling.checkText(name, r);
      }
    } catch (final IOException e) {
      Assert.fail(e.getMessage());
    }
  }

}

