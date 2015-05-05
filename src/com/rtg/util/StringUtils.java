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
package com.rtg.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableSet;

/**
 * String utilities.
 *
 */
public final class StringUtils {

  /** System dependent line separator. */
  public static final String LS = System.lineSeparator();
  /**  System dependent file separator. */
  public static final String FS = System.getProperty("file.separator");
  /** Tab. */
  public static final String TAB = "\t";

  private static final Map<Character, String> REPLACE = new HashMap<>();
  //see http://www.java-tips.org/java-se-tips/java.lang/character-escape-codes-in-java.html
  static {
    REPLACE.put('\n', "n");
    REPLACE.put('\t', "t");
    REPLACE.put('\b', "b");
    REPLACE.put('\r', "r");
    REPLACE.put('\f', "f");
    REPLACE.put('\\', "\\");
    REPLACE.put('\"', "\"");
    REPLACE.put('\'', "\'");
  }

  private StringUtils() { }

  /**
   * Converts line endings in the given string to line endings of current platform
   * @param input string to convert
   * @return converted string
   */
  public static String convertLineEndings(final String input) {
    //    final char[] inChar = input.toCharArray();
    //    final StringBuilder ret = new StringBuilder(input.length());
    //    char temp = '\0';
    //    int last = 0;
    //    for (int i = 0; i < inChar.length; i++) {
    //      if (inChar[i] == '\n' || inChar[i] == '\r') {
    //        if (temp == '\0' || temp == inChar[i]) {
    //          ret.append(new String(inChar, last, i - last));
    //          ret.append(LS);
    //          temp = inChar[i];
    //        } else {
    //          temp = '\0';
    //        }
    //        last = i + 1;
    //      }
    //    }
    //    if (last < inChar.length) {
    //      ret.append(inChar, last, inChar.length - last);
    //    }
    //    return ret.toString();
    //Is this too slow?
    return input.replaceAll("(\r\n)|\r|\n", LS);
  }

  /**
   * Return a string consisting of the specified number of spaces.
   *
   * @param length number of spaces
   * @return space string
   * @exception IllegalArgumentException if <code>length</code> is negative.
   */
  public static String getSpaceString(final int length) {
    return getCharString(' ', length);
  }

  /**
   * Return a string consisting of the specified character, repeated.
   *
   *
   * @param theChar the character to repeat
   * @param length number of repetitions
   * @return the string
   * @exception IllegalArgumentException if <code>length</code> is negative.
   */
  public static String getCharString(char theChar, final int length) {
    if (length < 0) {
      throw new IllegalArgumentException("bad length");
    }
    final char[] buf = new char[length];
    Arrays.fill(buf, theChar);
    return new String(buf);
  }

  /**
   * Faster implementation of string splitting on a single character delimiter.
   * @param src the string to split
   * @param delim the delimiter character
   * @return a <code>String[]</code> value
   */
  public static String[] split(String src, char delim) {
    return split(src, delim, 0);
  }

  /**
   * Faster implementation of string splitting on a single character delimiter.
   * @param src the string to split
   * @param delim the delimiter character
   * @param limit return at most limit - 1 entries, 0 or below for as many as necessary.
   * @return a <code>String[]</code> value
   */
  public static String[] split(String src, char delim, int limit) {
    final ArrayList<String> output = new ArrayList<>();
    int index;
    int lindex = 0;
    final boolean limited = limit > 0;
    while ((index = src.indexOf(delim, lindex)) != -1 && (!limited || output.size() < limit - 1)) {
      output.add(src.substring(lindex, index));
      lindex = index + 1;
    }
    output.add(src.substring(lindex));
    return output.toArray(new String[output.size()]);
  }

  /**
   * Format a number with commas.
   * @param n number to be formatted.
   * @return string with the formatted number.
   */
  public static String commas(final long n) {
    return String.format("%,d", n);
  }

  /**
   * Check if a string is a valid java identifier.
   * Spaces are not allowed fore or aft.
   * @param str to be tested.
   * @return true iff the string is a valid java identifier.
   */
  public static boolean isJavaIdentifier(final String str) {
    if (str == null || str.length() < 1) {
      return false;
    }
    if (!Character.isJavaIdentifierStart(str.charAt(0))) {
      return false;
    }
    for (int i = 1; i < str.length(); i++) {
      if (!Character.isJavaIdentifierPart(str.charAt(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Generate a string that contains <code>rep</code> repeated n times.
   * Returns length 0 string if n is less than 0.
   * @param rep string to be repeated.
   * @param n number of times to repeat it.
   * @return the repetitious string.
   */
  public static String repeat(final String rep, final int n) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; i++) {
      sb.append(rep);
    }
    return sb.toString();
  }

  /**
   * Translate a string into a form where hidden escape characters etc. are explicitly displayed.
   * The result should be able to be used as a literal in Java code. (it includes the quotes at the start and
   * end).
   * @param arg string to be displayed.
   * @return transformed string.
   */
  public static String display(final String arg) {
    final StringBuilder sb = new StringBuilder();
    sb.append("\"");
    for (int i = 0; i < arg.length(); i++) {
      final Character c = arg.charAt(i);
      final String trans = REPLACE.get(c);
      if (trans != null) {
        sb.append("\\").append(trans);
        continue;
      }
      if (c >= 32 && c < 127) {
        //in ASCII range
        sb.append(c);
        continue;
      }
      //use unicode escape for everything else
      //make sure 4 characters emitted
      sb.append("\\u");
      final String hexl = Integer.toHexString(c);
      final String hex = hexl.toUpperCase(Locale.getDefault());
      for (int j = hex.length(); j < 4; j++) {
        sb.append("0");
      }
      sb.append(hex);
    }
    sb.append("\"");
    return sb.toString();
  }

  /**
   * Pad a string on the left with spaces till it is at least the specified length.
   * @param s initial string.
   * @param length to pad the string to.
   * @return padded string.
   */
  public static String padLeft(final String s, final int length) {
    final int pad = length - s.length();
    return pad > 0 ? getSpaceString(pad) + s : s;
  }

  /**
   * Pad a string on the right with spaces till it is at least the specified length.
   * @param s initial string.
   * @param length to pad the string to.
   * @return padded string.
   */
  public static String padRight(final String s, final int length) {
    final int pad = length - s.length();
    return pad > 0 ? s + getSpaceString(pad): s;
  }

  /**
   * Pads out two string by inserting spaces between them.
   *
   * @param first left string
   * @param length total length of the desired output string
   * @param last right string
   * @return the padded string
   */
  public static String padBetween(final String first, final int length, final String last) {
    final int padlen = length - (first.length() + last.length());
    return first + (padlen > 0 ? getSpaceString(padlen) : "") + last;
  }


  /**
   * Implodes a string array into a String.
   * @param array the array to implode
   * @param delim the delimiter to insert between values
   * @param quoteStringsWithSpaces whether to surround strings that contain spaces with quotes or not
   * @return the imploded string
   */
  public static String implode(final String[] array, final String delim, final boolean quoteStringsWithSpaces) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < array.length; i++) {
      if (i != 0 && delim != null) {
        sb.append(delim);
      }
      sb.append(quoteStringsWithSpaces ? smartQuote(array[i]) : array[i]);
    }
    return sb.toString();
  }

  /**
   * Implodes a string array into a String.
   * @param array the array to implode
   * @param delim the delimiter to insert between values
   * @return the imploded string
   */
  public static String implode(final String[] array, final String delim) {
    return implode(array, delim, false);
  }

  /**
   * Implodes a collection of strings into a String.
   * @param strings the collection to implode
   * @param delim the delimiter to insert between values
   * @param quoteStringsWithSpaces whether to surround strings that contain spaces with quotes or not
   * @return the imploded string
   */
  public static String implode(final Collection<String> strings, final String delim, final boolean quoteStringsWithSpaces) {
    return implode(strings.toArray(new String[strings.size()]), delim, quoteStringsWithSpaces);
  }

  /**
   * Implodes a collection of strings into a String.
   * @param strings the collection to implode
   * @param delim the delimiter to insert between values
   * @return the imploded string
   */
  public static String implode(final Collection<String> strings, final String delim) {
    return implode(strings, delim, false);
  }

  /**
   * Converts a two's complement hex string to long,
   * i.e. reverses {@link java.lang.Long#toHexString(long)}
   * @param hex two's complement hex string
   * @return long value
   */
  public static long fromLongHexString(String hex) {
    long ret = 0;
    for (final char c : hex.toCharArray()) {
      final int dig = Character.digit(c, 16);
      if (dig == -1) {
        throw new NumberFormatException(c + " is not a valid hexadecimal digit");
      }
      ret = (ret << 4) | dig;
    }
    return ret;
  }


  /**
   * Makes the given string safe for use as XML character data. This
   * involves replacing <code>&lt;</code>, <code>&gt;</code>,
   * <code>&amp;</code> and <code>&quot;</code> with protected equivalents
   * <code>&amp;lt;</code>, <code>&amp;gt;</code> and <code>&amp;amp;</code>
   * and <code>&amp;quot;</code>
   *
   * @param s string to protect
   * @return protected version of s
   */
  public static String xmlProtect(final String s) {
    if (s == null) {
      return null;
    }
    final int slen = s.length();
    if (slen == 0) {
      return s;
    }
    final StringBuilder b = new StringBuilder();
    for (int i = 0; i < slen; i++) {
      final char c = s.charAt(i);
      switch (c) {
        case '<':
          b.append("&lt;");
          break;
        case '>':
          b.append("&gt;");
          break;
        case '&':
          b.append("&amp;");
          break;
        case '"':
          b.append("&quot;");
          break;
        default:
          b.append(c);
          break;
      }
    }
    return b.toString();
  }

  /**
   * Like grep, but in java. Unfortunately this will replace any new lines with the current platforms new line.
   * @param val value to search
   * @param pattern pattern to match
   * @return lines containing pattern
   */
  public static String grep(String val, String pattern) {
    return grepInner(val, pattern, false);
  }
  private static String grepInner(String val, String pattern, boolean minusV) {
    final StringBuilder ret = new StringBuilder();
    try {
      try (BufferedReader r = new BufferedReader(new StringReader(val))) {
        String line;
        while ((line = r.readLine()) != null) {
          if (line.matches(".*(" + pattern + ").*") ^ minusV) {
            ret.append(line).append(LS);
          }
        }
      }
    } catch (final IOException e) {
      //pretty sure this is impossible without programmer error current setup (reading from string source)
      throw new IllegalStateException("Improper use of StringReader", e);
    }
    return ret.toString();
  }
  /**
   * Like grep -v, but in java. Unfortunately this will replace any new lines with the current platforms new line.
   * @param val value to search
   * @param pattern pattern to match
   * @return lines not containing pattern
   */
  public static String grepMinusV(String val, String pattern) {
    return grepInner(val, pattern, true);
  }

  /**
   * Produces a string that consists of the <code>toString</code> of each element in collection with the glue string separating each item
   * @param <T> type of object to join
   * @param glue the string that will be inserted between elements of items
   * @param items the elements to join into a string
   * @return a string containing all elements of items separated by glue
   */
  public static <T> String join(String glue, Collection<T> items) {
    final StringBuilder sb = new StringBuilder();
    String join = "";
    for (final T item : items) {
      sb.append(join);
      join = glue;
      sb.append(item.toString());
    }
    return sb.toString();
  }


  private static boolean equalsLeft(final String a, final String b, int pos, final int rightOffset) {
    return pos < (a.length() - rightOffset) && pos < (b.length() - rightOffset) && a.charAt(pos) == b.charAt(pos);
  }

  private static boolean equalsRight(final String a, final String b, int pos, final int leftOffset) {
    final int lpos = pos + leftOffset;
    return lpos < a.length() && lpos < b.length() && a.charAt(a.length() - pos - 1) == b.charAt(b.length() - pos - 1);
  }

  /**
   * Return the length of the longest common prefix of the supplied strings.
   * @param strings strings to test
   * @return longest common prefix
   */
  public static int longestPrefix(final String[] strings) {
    return longestPrefix(strings, 0);
  }
  /**
   * Return the length of the longest common prefix of the supplied strings.
   * @param strings strings to test
   * @param rightOffset effective right edge of strings (i.e. do not find prefix going into this region)
   * @return longest common prefix
   */
  public static int longestPrefix(final String[] strings, final int rightOffset) {
    if (strings.length <= 1) {
      return strings.length == 0 ? 0 : strings[0].length();
    }
    final String a = strings[0];
    int clip = -1;
    while (true) {
      clip++;
      for (int k = 1; k < strings.length; k++) {
        if (!equalsLeft(a, strings[k], clip, rightOffset)) {
          return clip;
        }
      }
    }
  }

  /**
   * Return the length of the longest common suffix of the supplied strings.
   * @param strings strings to test
   * @param leftOffset effective left edge of strings (i.e. do not find suffix going into this region)
   * @return longest common suffix
   */
  public static int longestSuffix(final String[] strings, final int leftOffset) {
    if (strings.length <= 1) {
      return strings.length == 0 ? 0 : strings[0].length();
    }
    final String a = strings[0];
    int clip = -1;
    while (true) {
      clip++;
      for (int k = 1; k < strings.length; k++) {
        if (!equalsRight(a, strings[k], clip, leftOffset)) {
          return clip;
        }
      }
    }
  }

  /**
   * Return the length of the longest common suffix of the supplied strings.
   * @param strings strings to test
   * @return longest common suffix
   */
  public static int longestSuffix(final String[] strings) {
    return longestSuffix(strings, 0);
  }

  /**
   * Remove backslash escape sequences from a string.
   * @param s string to remove escapes from
   * @return string without backslash escapes
   */
  public static String removeBackslashEscapes(final String s) {
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < s.length(); k++) {
      final char c = s.charAt(k);
      if (c == '\\') {
        if (++k == s.length()) {
          break;
        }
        final char d = s.charAt(k);
        switch (d) {
          case 'n':
            sb.append('\n');
            break;
          case 'r':
            sb.append('\r');
            break;
          case 't':
            sb.append('\t');
            break;
          default: // \, #, etc.
            sb.append(d);
            break;
        }
      } else if (c == '#') {
        break;
      } else {
        sb.append(c);
      }
    }
    return sb.toString();
  }

  /**
   * Does an deep copy of the input string and returns so <code>input != (return value)</code>
   * @param s input string
   * @return deep copy
   */
  public static String deepCopy(final String s) {
    return new String(s);
  }

  /**
   * Quotes a string, escaping embedded slash or quotes
   * @param s string to quote
   * @return quoted string, always
   */
  public static String dumbQuote(String s) {
    return "\"" + s.replaceAll("\\\\", "\\\\\\\\").replaceAll("\"", "\\\\\"") + "\"";
  }

  /**
   * Quotes string when spaces are detected
   * @param s string to quote
   * @return quoted string, or not
   */
  public static String smartQuote(String s) {
    if (s.contains(" ")) {
      return dumbQuote(s);
    } else {
      return s;
    }
  }

  /**
   * Expand a name if it is an unambiguous prefix of an entry in the given set of names
   * @param names the reference set of names
   * @param prefix the name prefix
   * @return the input name, or an unambiguous expansion of it
   */
  public static String expandPrefix(NavigableSet<String> names, String prefix) {
    String name = prefix;
    if (!names.contains(name)) {
      // Try a search for unambiguous prefix
      final String closest = names.ceiling(prefix);
      if (closest != null && closest.startsWith(prefix)) {
        final String higher = names.higher(closest);
        if (higher == null || !higher.startsWith(prefix)) {
          name = closest;
        }
      }
    }
    return name;
  }

  /**
   * Reverse the string.
   * @param str to be reversed.
   * @return reversed string.
   */
  public static String reverse(final String str) {
    if (str == null) {
      return null;
    }
    return new StringBuilder(str).reverse().toString();
  }
}
