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

package com.rtg.visualization;



/**
 * This class is concerned with display of text with <code>ANSI</code> color markup
 */
public class HtmlDisplayHelper extends DisplayHelper {

  static final char TAG_START = '<';
  static final char TAG_END = '>';

  static final String END_SPAN = "</span>";
  static final String BOLD_ON = "<span style=\"font-weight:bold\">";
  static final String UNDERLINE_ON = "<span style=\"text-decoration:underline\">";

  static final String[] FGCOLORS = {
    "<span style=\"color:black\">",
    "<span style=\"color:red\">",
    "<span style=\"color:green\">",
    "<span style=\"color:orange\">",
    "<span style=\"color:blue\">",
    "<span style=\"color:magenta\">",
    "<span style=\"color:darkcyan\">",
    "<span style=\"color:white\">",
    "<span style=\"color:whitesmoke\">"
  };

  static final String[] BGCOLORS = {
    "<span style=\"background-color:black\"",
    "<span style=\"background-color:red\">",
    "<span style=\"background-color:green\">",
    "<span style=\"background-color:orange\">",
    "<span style=\"background-color:dodgerblue\">",
    "<span style=\"background-color:magenta\">",
    "<span style=\"background-color:cyan\">",
    "<span style=\"background-color:white\">",
    "<span style=\"background-color:whitesmoke\">"
  };

  @Override
  protected String header() {
    return "<html><body><pre>";
  }
  @Override
  protected String footer() {
    return "</pre></body></html>";
  }

  @Override
  public boolean supportsNesting() {
    return true;
  }

  @Override
  public String decorateUnderline(String text) {
    return UNDERLINE_ON + text + END_SPAN;
  }

  @Override
  public String decorateBold(String text) {
    return BOLD_ON + text + END_SPAN;
  }

  @Override
  public String decorateForeground(String text, int color) {
    return FGCOLORS[color] + text + END_SPAN;
  }

  @Override
  public String decorateBackground(String text, int color) {
    return BGCOLORS[color] + text + END_SPAN;
  }

  @Override
  protected boolean isMarkupStart(char c) {
    return c == TAG_START;
  }
  @Override
  protected boolean isMarkupEnd(char c) {
    return c == TAG_END;
  }

  @Override
  protected String escape(String text) {
    final StringBuilder sb = new StringBuilder();
    for (int currpos = 0; currpos < text.length(); currpos++) {
      final char c = text.charAt(currpos);
      if (c == TAG_START) {
        sb.append("&lt;");
      } else if (c == TAG_END) {
        sb.append("&gt;");
      } else {
        sb.append(c);
      }
    }
    return sb.toString();
  }
}
