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
public class AnsiDisplayHelper extends DisplayHelper {

  static final char ESC = (char) 27;
  static final int DEFAULT = 9;

  static final String RESET_ATTRIBUTES = ESC + "[0m";
  static final String BOLD_ON = ESC + "[1m";
  static final String BOLD_OFF = ESC + "[22m";
  static final String DIM_ON = ESC + "[2m"; /// Not widely supported
  static final String DIM_OFF = ESC + "[22m";
  static final String ITALICS_ON = ESC + "[3m"; /// Not widely supported, sometimes shown as inverse
  static final String ITALICS_OFF = ESC + "[23m"; /// Not widely supported, sometimes shown as inverse
  static final String UNDERLINE_ON = ESC + "[4m";
  static final String UNDERLINE_OFF = ESC + "[24m";
  static final String INVERSE_ON = ESC + "[7m";
  static final String INVERSE_OFF = ESC + "[27m";

  static int ansiColor(int color) {
    if (color == WHITE_PLUS) {
      return DEFAULT;
    }
    return color;
  }

  static String ansiForeground(int color) {
    return ESC + "[3" + ansiColor(color) + "m";
  }

  static String ansiBackground(int color) {
    return ESC + "[4" + ansiColor(color) + "m";
  }

  @Override
  public String decorateUnderline(String text) {
    return UNDERLINE_ON + text + UNDERLINE_OFF;
  }

  @Override
  public String decorateBold(String text) {
    return BOLD_ON + text + BOLD_OFF;
  }

  @Override
  public String decorateForeground(String text, int color) {
    return ansiForeground(color) + text + ansiForeground(DEFAULT);
  }

  @Override
  public String decorateBackground(String text, int color) {
    return ansiBackground(color) + text + ansiBackground(DEFAULT);
  }

  @Override
  protected boolean isMarkupStart(char c) {
    return c == ESC;
  }
  @Override
  protected boolean isMarkupEnd(char c) {
    return c == 'm';
  }
}
