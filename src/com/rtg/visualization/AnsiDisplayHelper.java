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

  static final String[] FGCOLORS = {
    ESC + "[3" + BLACK + "m",
    ESC + "[3" + RED + "m",
    ESC + "[3" + GREEN + "m",
    ESC + "[3" + YELLOW + "m",
    ESC + "[3" + BLUE + "m",
    ESC + "[3" + MAGENTA + "m",
    ESC + "[3" + CYAN + "m",
    ESC + "[3" + WHITE + "m",
    ESC + "[39m"
  };

  // These background colors require 256 color support
  static final String[] BGCOLORS = {
    ESC + "[48;5;" + extendedColor(0, 0, 0) + "m",
    ESC + "[48;5;" + extendedColor(1, 0, 0) + "m",
    ESC + "[48;5;" + extendedColor(0, 1, 0) + "m",
    ESC + "[48;5;" + extendedColor(1, 1, 0) + "m",
    ESC + "[48;5;" + extendedColor(0, 0, 1) + "m",
    ESC + "[48;5;" + extendedColor(1, 0, 1) + "m",
    ESC + "[48;5;" + extendedColor(0, 1, 1) + "m",
    ESC + "[48;5;" + extendedColor(1, 1, 1) + "m",
    ESC + "[48;5;237m",
    //ESC + "[49m"
  };


  // Map from 0-6 RGB values into ANSI 256 color extended range
  static int extendedColor(int r, int g, int b) {
    if ((r >= 6) || (g >= 6) || (b >= 6)) {
      throw new IllegalArgumentException();
    }
    return 16 + r * 36 + g * 6 + b;
  }

  static String defaultForeground() {
    return ESC + "[39m";
  }
  static String defaultBackground() {
    return ESC + "[49m";
  }

  static String ansiForeground(int color) {
    return FGCOLORS[color];
  }

  static String ansiBackground(int color) {
    return BGCOLORS[color];
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
    return ansiForeground(color) + text + defaultForeground();
  }

  @Override
  public String decorateBackground(String text, int color) {
    return ansiBackground(color) + text + defaultBackground();
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
