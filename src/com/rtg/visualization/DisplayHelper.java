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


import com.rtg.util.StringUtils;

/**
 * This class is concerned with display of text without any markup
 */
public class DisplayHelper {

  static final int LABEL_LENGTH = 6;

  static final char SPACE_CHAR = ' ';
  static final char INSERT_CHAR = '_';

  //Define color codes. These currently must match ANSI color codes, see http://en.wikipedia.org/wiki/ANSI_escape_code#CSI_Codes
  /** Color code for black */
  public static final int BLACK = 0;
  /** Color code for red */
  public static final int RED = 1;
  /** Color code for green */
  public static final int GREEN = 2;
  /** Color code for yellow */
  public static final int YELLOW = 3;
  /** Color code for blue */
  public static final int BLUE = 4;
  /** Color code for magenta */
  public static final int MAGENTA = 5;
  /** Color code for cyan */
  public static final int CYAN = 6;
  /** Color code for white */
  public static final int WHITE = 7;

  /** Color code for off-white */
  public static final int WHITE_PLUS = 8;


  String getSpaces(final int diff) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < diff; i++) {
      sb.append(SPACE_CHAR);
    }
    return sb.toString();
  }

  String getInserts(final int diff) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < diff; i++) {
      sb.append(INSERT_CHAR);
    }
    return sb.toString();
  }

  int getBaseColor(char readChar) {
    switch (readChar) {
      case 'a':
      case 'A':
        return DisplayHelper.GREEN;
      case 't':
      case 'T':
        return DisplayHelper.RED;
      case 'c':
      case 'C':
        return DisplayHelper.BLUE;
      case 'g':
      case 'G':
        return DisplayHelper.MAGENTA;
      default:
        return -1;
    }
  }

  protected String escape(String text) {
    return text;
  }
  protected String header() {
    return null;
  }
  protected String footer() {
    return null;
  }

  /**
   * @return true if display properties can nest
   */
  public boolean supportsNesting() {
    return false;
  }

  /**
   * Mark up the text with underlining
   * @param text the text
   * @return the marked up text
   */
  public String decorateUnderline(String text) {
    return text;
  }

  /**
   * Mark up the text with boldness
   * @param text the text
   * @return the marked up text
   */
  public String decorateBold(String text) {
    return text;
  }

  /**
   * Mark up the text with a foreground color
   * @param text the text
   * @param color the foreground color
   * @return the marked up text
   */
  public String decorateForeground(String text, int color) {
    return text;
  }

  /**
   * Mark up the text with a background color
   * @param text the text
   * @param color the background color
   * @return the marked up text
   */
  public String decorateBackground(String text, int color) {
    return text;
  }

  /**
   * Mark up the text with foreground and background colors
   * @param text the text
   * @param fgcolor the foreground color
   * @param bgcolor the background color
   * @return the marked up text
   */
  public String decorate(final String text, int fgcolor, int bgcolor) {
    return decorateForeground(decorateBackground(text, bgcolor), fgcolor);
  }

  protected String decorateLabel(final String label) {
    String shortLabel = (label.length() >= LABEL_LENGTH ? label.substring(label.length() - LABEL_LENGTH) : label) + ":";
    if (shortLabel.length() < LABEL_LENGTH) {
      shortLabel = getSpaces(LABEL_LENGTH - shortLabel.length()) + shortLabel;
    }
    return decorateForeground(shortLabel, DisplayHelper.CYAN) + " ";
  }

  // Decorates a section of DNA with highlight colors, not marking up space characters at the ends
  protected String trimHighlighting(final String dna, int bgcolor) {
    final String trimmed = StringUtils.trimSpaces(dna);
    if (trimmed.length() == 0) { // All whitespace, no markup needed
      return dna;
    }
    if (trimmed.length() == dna.length()) { // No trimming needed
      return decorateBackground(dna, bgcolor);
    } else {
      if (dna.charAt(0) == ' ') { // Some amount of non-marked up prefix needed
        int prefixEnd = 0;
        while (dna.charAt(prefixEnd) == ' ') {
          prefixEnd++;
        }
        return dna.substring(0, prefixEnd) + decorateBackground(trimmed, bgcolor) + dna.substring(prefixEnd + trimmed.length());
      } else { // Trimming was only at the end
        return decorateBackground(trimmed, bgcolor) + dna.substring(trimmed.length());
      }
    }
  }

  protected boolean isMarkupStart(char c) {
    return false;
  }
  protected boolean isMarkupEnd(char c) {
    return false;
  }

  protected String decorateWithHighlight(final String str, boolean[] highlightMask, int bgcolor, boolean colorBases) {
    final StringBuilder output = new StringBuilder();
    int coord = 0; // coordinate ignoring markup
    final StringBuilder toHighlight = new StringBuilder();
    boolean highlight = false;
    boolean inMarkup = false;
    for (int i = 0; i < str.length(); i++) {
      final char c = str.charAt(i);
      if (inMarkup) {
        if (isMarkupEnd(c)) {
          inMarkup = false;
        }
        output.append(c);
      } else {
        if (isMarkupStart(c)) {
          inMarkup = true;
          output.append(c);
        } else {
          if (highlightMask != null && highlightMask[coord] != highlight) {
            if (highlight) {
              output.append(trimHighlighting(toHighlight.toString(), bgcolor));
              toHighlight.setLength(0);
            }
            highlight = highlightMask[coord];
          }
          final StringBuilder dest = highlight ? toHighlight : output;
          final int col = getBaseColor(c);
          if (colorBases && col >= 0) {
            dest.append(decorateForeground(String.valueOf(c), col));
          } else {
            dest.append(c);
          }
          coord++;
        }
      }
    }
    if (highlight) {
      output.append(trimHighlighting(toHighlight.toString(), bgcolor));
    }
    return output.toString();
  }

  /**
   * Trims a sequence for display.
   * @param sequence the string to clip. May contain markup
   * @param clipStart first position in non-markup coordinates to output
   * @param clipEnd end position (exclusive) in non-markup coordinates
   * @return the clipped sequence.
   */
  protected String clipSequence(String sequence, final int clipStart, final int clipEnd) {
    final StringBuilder sb = new StringBuilder();
    int coord = 0; // coordinate ignoring markup
    boolean inMarkup = false;
    for (int currpos = 0; currpos < sequence.length(); currpos++) {
      final char c = sequence.charAt(currpos);
      if (inMarkup) {
        if (isMarkupEnd(c)) {
          inMarkup = false;
        }
        sb.append(c);
      } else {
        if (isMarkupStart(c)) {
          inMarkup = true;
          sb.append(c);
        } else {
          if ((coord >= clipStart) && (coord < clipEnd)) {
            sb.append(c);
          }
          coord++;
        }
      }
    }
    return sb.toString();
  }
}
