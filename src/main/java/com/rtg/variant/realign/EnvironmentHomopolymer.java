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
    for (int i = -maxShift(); i < readLength() + maxShift(); ++i) {
      final int j = absoluteTemplatePosition(i);
      if (j >= 0 && j < templateLength()) {
        sb.append("[").append(i).append("]").append(VALUE_CHARS[template(i)]).append(" ").append(templateStart(i)).append(" ").append(templateEnd(i)).append(LS);
      }
    }
    sb.append("Read 0..").append(readLength()).append(LS);
    //System.err.println("readLength()=" + readLength());
    for (int i = 0; i < readLength(); ++i) {
      sb.append("[").append(i).append("]").append(VALUE_CHARS[read(i)]).append(" ").append(Utils.realFormat(quality(i), 3)).append(" ").append(readStart(i)).append(" ").append(readEnd(i)).append(LS);
    }
    return sb.toString();
  }
}
