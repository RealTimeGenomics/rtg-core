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

package com.rtg.reference;

import static com.rtg.util.StringUtils.LS;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Locale;
import java.util.Map;
import java.util.Queue;

import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Parsing for reference file.
 */
class ReferenceParse {

  static RegionRestriction region(final Map<String, Integer> sequenceLengths, final String str) {
    try {
      final RegionRestriction reg = new RegionRestriction(str);
      final String name = reg.getSequenceName();
      final Integer length = sequenceLengths.get(name);

      if (length == null || reg.getEnd() > length) {
        return null;
      }
      return reg;
    } catch (final IllegalArgumentException e) {
      return null;
    }
  }

  static Ploidy getPloidy(final String ploidyString) {
    try {
      return Ploidy.valueOf(ploidyString.toUpperCase(Locale.getDefault()));
    } catch (final IllegalArgumentException e) {
      return null;
    }
  }

  static Sex getSex(final String sexStr) {
    try {
      return Sex.valueOf(sexStr.toUpperCase(Locale.getDefault()));
    } catch (final IllegalArgumentException e) {
      return null;
    }
  }

  /**
   * Test if lengths of regions the same.
   * @param r1 first region to be checked.
   * @param r2 second region to be tested.
   * @return true iff their lengths are the same.
   */
  static boolean lengthsEqual(final Map<String, Integer> sequenceLengths, RegionRestriction r1, RegionRestriction r2) {
    return length(sequenceLengths, r1) == length(sequenceLengths, r2);
  }

  static int length(final Map<String, Integer> names, RegionRestriction r) {
    final int start = r.getStart() == RegionRestriction.MISSING ? 0 : r.getStart();
    final int end = r.getEnd() == RegionRestriction.MISSING ? names.get(r.getSequenceName()) : r.getEnd();
    return end - start;
  }

  static boolean sexMatch(final Sex actual, final Sex lineSex) {
    return lineSex == Sex.EITHER || actual == lineSex;
  }

  static Boolean linear(final String linearString) {
    switch (linearString) {
      case "linear":
        return true;
      case "circular":
        return false;
      default:
        return null;
    }
  }

  /**
   * Check for comments, trim white space. Return null if nothing left after that.
   * Otherwise split on tabs and return the array.
   * @param line to be split.
   * @return the split line or null if a problem.
   */
  static String[] splitLine(final String line) {
    final int ix0 = line.indexOf('#');
    final int ix = ix0 == -1 ? line.length() : ix0;
    final String lessComment = line.substring(0, ix);
    if (lessComment.matches("^\\s*$")) {
      return null;
    }
    return lessComment.split("\\s+");
  }

  final Map<String, ReferenceSequence> mReferences = new HashMap<>();
  final Map<String, Integer> mNames;
  final BufferedReader mIn;
  final Sex mSex;
  boolean mError = false;
  final Queue<Pair<RegionRestriction, RegionRestriction>> mDuplicates = new LinkedList<>();
  boolean mLinearDefault;
  Ploidy mPloidyDefault = null;
  int mNonblankLines = 0;

  ReferenceParse(Map<String, Integer> names, BufferedReader in, Sex sex) {
    mNames = names;
    mIn = in;
    this.mSex = sex;
  }

  void error(final String msg) {
    Diagnostic.warning(msg);
    mError = true;
  }

  /**
   * Parse a reference file putting results into fields.
   * @throws IOException if I/O exception reading reference file.
   */
  void parse() throws IOException {
    while (true) {
      final String line = mIn.readLine();
      if (line == null) {
        break;
      }
      final String msg = line(line);
      if (msg != null) {
        error("Error reading reference file on line:" + line + LS + msg);
      }
    }
    end();
  }

  /**
   * Parse one line of reference file.
   * @param line complete non-null line.
   * @return an error message or null if no error.
   */
  String line(final String line) {

    final String[] split = splitLine(line);
    if (split == null) {
      return null;
    }
    mNonblankLines++;
    if (split.length < 2) {
      return "Version line too short";
    }
    if (mNonblankLines == 1) {
      //version line
      if (split.length != 2 || !split[0].equals("version")
          || !(split[1].equals("0") || split[1].equals("1"))) {
        return "Invalid version line.";
      }
      return null;
    }

    //normal lines
    final Sex lineSex = getSex(split[0]);
    if (lineSex == null) {
      return "Invalid sex:" + split[0];
    }
    final String type = split[1];

    final boolean match = sexMatch(mSex, lineSex);
    switch (type) {
      case "def":
        return def(split, match);
      case "seq":
        return seq(split, match);
      case "dup":
        return dup(split, match);
      default:
        return "Invalid line type (should be one of: def, seq, dup):" + type;
    }
  }

  /**
   * Parse a duplicate line.
   * @param split line after splitting on tabs.
   * @param match if true sexes match and actually carry out operation.
   * @return an error message or null if no errors.
   */
  String dup(final String[] split, final boolean match) {
    if (split.length != 4) {
      return "Duplicate line has incorrect number of fields.";
    }
    final String dup1 = split[2];
    final RegionRestriction r1 = region(mNames, dup1);
    if (r1 == null) {
      return "Invalid region:" + dup1;
    }
    final String dup2 = split[3];
    final RegionRestriction r2 = region(mNames, dup2);
    if (r2 == null) {
      return "Invalid region:" + dup2;
    }
    if (!lengthsEqual(mNames, r1, r2)) {
      return "Lengths of regions disagree.";
    }
    if (match) {
      mDuplicates.add(new Pair<>(r1, r2));
    }
    return null;
  }

  /**
   * Parse a sequence line.
   * @param split line after splitting on tabs.
   * @param match if true sexes match and actually carry out operation.
   * @return an error message or null if no errors.
   */
  String seq(final String[] split, final boolean match) {
    if (split.length < 5) {
      return "Sequence line has incorrect number of fields.";
    }
    final String name = split[2].trim();
    if ("".equals(name) || name.contains(" ")) {
      return "Invalid sequence name:" + name;
    }
    if (!mNames.containsKey(name)) {
      return "Sequence in reference file:" + name + " not found in genome.";
    }
    final String ploidyString = split[3];
    final Ploidy ploidy = getPloidy(ploidyString);
    if (ploidy == null) {
      return "Invalid ploidy value:" + ploidyString;
    }
    final String linearString = split[4];
    final Boolean linear = linear(linearString);
    if (linear == null) {
      return "Invalid linear/circular value:" + linearString;
    }
    final String hapMate; // haploid complement is optional e.g. X0 sex determination system.
    if ((ploidy == Ploidy.HAPLOID) && (split.length == 6)) {
      hapMate = split[5].trim();
      if ("".equals(hapMate) || hapMate.contains(" ")) {
        return "Invalid haploid mate sequence name:" + hapMate;
      }
      if (!mNames.containsKey(hapMate)) {
        return "Haploid mate sequence in reference file:" + hapMate + " not found in genome.";
      }
    } else if (split.length > 5) {
      return "Sequence line has incorrect number of fields.";
    } else {
      hapMate = null;
    }

    if (match) {
      if (mReferences.containsKey(name)) {
        return "Sequence defined twice:" + name;
      }
      mReferences.put(name, new ReferenceSequence(true, linear, ploidy, name, hapMate, mNames.get(name)));
    }
    return null;
  }

  /**
   * Parse a default line.
   * @param split line after splitting on tabs.
   * @param match if true sexes match and actually carry out operation.
   * @return an error message or null if no errors.
   */
  String def(final String[] split, final boolean match) {
    if (split.length != 4) {
      return "Default line has incorrect number of fields.";
    }
    final String ploidyString = split[2];
    final Ploidy ploidy = getPloidy(ploidyString);
    if (ploidy == null) {
      return "Invalid ploidy value:" + ploidyString;
    }
    final String linearString = split[3];
    final Boolean linear = linear(linearString);
    if (linear == null) {
      return "Invalid linear/circular value:" + linearString;
    }
    if (match) {
      if (mPloidyDefault != null) {
        return "Duplicate default definition.";
      }
      mLinearDefault = linear;
      mPloidyDefault = ploidy;
    }
    return null;
  }

  public void end() {
    if (mNonblankLines == 0) {
      error("No valid lines found in reference file.");
    }
  }
}
