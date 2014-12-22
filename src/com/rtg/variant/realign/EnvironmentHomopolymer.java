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

package com.rtg.variant.realign;

import static com.rtg.util.StringUtils.LS;

import com.rtg.mode.DNA;
import com.rtg.util.Utils;


/**
 */
public class EnvironmentHomopolymer extends EnvironmentDecorator {

  private static final char[] VALUE_CHARS = DNA.valueChars();

  final HomopolymerRepeats mTemplate;

  final HomopolymerRepeats mRead;

  /**
   * @param env original environment.
   */
  public EnvironmentHomopolymer(final Environment env) {
    super(env);

    final ByteArrayAdaptor templateSequence = new ByteArrayAdaptor() {
      @Override
      public int length() {
        return mEnv.readLength() + 2 * mEnv.maxShift(); //yes it really is the read length
      }
      @Override
      public byte get(int index) {
        return mEnv.template(index - mEnv.maxShift());
      }
    };
    //System.err.println(templateSequence);
    mTemplate = new HomopolymerRepeats(templateSequence);
    final ByteArrayAdaptor readSequence = new ByteArrayAdaptor() {
      @Override
      public int length() {
        return mEnv.readLength();
      }
      @Override
      public byte get(int index) {
        return mEnv.read(index);
      }
    };
    //System.err.println(readSequence);
    mRead = new HomopolymerRepeats(readSequence);
  }

  int templateStart(final int i) {
    return mTemplate.reverse(i + mEnv.maxShift());
  }

  int templateEnd(final int i) {
    return mTemplate.forward(i + mEnv.maxShift());
  }

  int readStart(final int i) {
    return mRead.reverse(i);
  }

  int readEnd(final int i) {
    return mRead.forward(i);
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Environment").append(LS);
    sb.append("Template ").append(absoluteTemplatePosition(0)).append("..").append(absoluteTemplatePosition(readLength() - 1)).append(LS);
    for (int i = -maxShift(); i < readLength() + maxShift(); i++) {
      final int j = absoluteTemplatePosition(i);
      if (j >= 0 && j < templateLength()) {
        sb.append("[").append(i).append("]").append(VALUE_CHARS[template(i)]).append(" ").append(templateStart(i)).append(" ").append(templateEnd(i)).append(LS);
      }
    }
    sb.append("Read 0..").append(readLength()).append(LS);
    //System.err.println("readLength()=" + readLength());
    for (int i = 0; i < readLength(); i++) {
      sb.append("[").append(i).append("]").append(VALUE_CHARS[read(i)]).append(" ").append(Utils.realFormat(quality(i), 3)).append(" ").append(readStart(i)).append(" ").append(readEnd(i)).append(LS);
    }
    return sb.toString();
  }
}
