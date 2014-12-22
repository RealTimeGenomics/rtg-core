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

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;


/**
 * Encapsulates a single flag.
 */
@TestClass(value = {"com.rtg.util.cli.CFlagsTest"})
public class Flag implements Comparable<Flag> {

  private final Character mFlagChar;

  private final String mFlagName;

  private final String mFlagDescription;

  /** The maximum number of times the flag can occur. */
  private int mMaxCount;

  /** The minimum number of times the flag can occur. */
  private int mMinCount;

  private final Class<?> mParameterType;

  private final String mParameterDescription;

  private Object mParameterDefault = null;

  private String mCategory = null;

  private String mPsuedoMinMaxString = null;

  /** Optional list of valid values for the parameter. */
  private List<String> mParameterRange;

  private boolean mRangeList = false;

  /** Values supplied by the user */
  private List<Object> mParameter = new ArrayList<>();

  /**
   * Creates a new <code>Flag</code> for which the name must be supplied on
   * the command line.
   *
   * @param flagChar a <code>Character</code> which can be supplied by the
   * user as an abbreviation for flagName. May be null.
   * @param flagName a <code>String</code> which is the name that the user
   * specifies on the command line to denote the flag.
   * @param flagDescription a name used when printing help messages.
   * @param minCount the minimum number of times the flag must be specified.
   * @param maxCount the maximum number of times the flag can be specified.
   * @param paramType a <code>Class</code> denoting the type of values to be
   * accepted. Maybe null for "switch" type flags.
   * @param paramDescription a description of the meaning of the flag.
   * @param paramDefault a default value that can be used for optional flags.
   * @param category The flag category
   */
  public Flag(final Character flagChar, final String flagName, final String flagDescription,
      final int minCount, final int maxCount, final Class<?> paramType, final String paramDescription,
      final Object paramDefault, final String category) {
    if (flagDescription == null) {
      throw new NullPointerException();
    }
    if (flagName == null) {
      if (paramType == null) {
        throw new IllegalArgumentException();
      }
    } else {
      if (flagName.startsWith("-")) {
        throw new IllegalArgumentException("Long flag names cannot start with '-'");
      }
    }
    setCategory(category);
    mFlagName = flagName;
    mFlagChar = flagChar;
    mFlagDescription = flagDescription;

    mParameterType = paramType;
    mParameterDescription = (mParameterType == null) ? null
        : ((paramDescription == null) || (paramDescription.length() == 0)) ? autoDescription(mParameterType)
            : paramDescription.toUpperCase(Locale.getDefault());
    if (mParameterType != null) {
      setParameterDefault(paramDefault);
    }

    mMinCount = minCount;
    mMaxCount = maxCount;

    // For enums set up the limited set of values message
    final String[] range = values(mParameterType);
    if (range != null) {
      setParameterRange(range);
    }
  }

  /**
   * Gets the object specified by str.
   * @param type of object.
   * @param str string to specify object.
   * @return object of class type (null if the type does not look sufficiently like an Enum).
   */
  static Object valueOf(final Class<?> type, final String str) {
    if (!isValidEnum(type)) {
      return null;
    }
    try {
      final String valueOfMethod = "valueOf";
      final Method m = type.getMethod(valueOfMethod, String.class);
      if (!isStatic(m)) {
        return null;
      }
      final Class<?> returnType = m.getReturnType();
      if (!type.isAssignableFrom(returnType)) {
        return null;
      }
      return m.invoke(null, str);
    } catch (final NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
      //Should never happen
      throw new RuntimeException(e);
    }
  }

  /**
   * Gets the range of values for an Enum (or at least something that looks like an Enum).
   * @param type from which values to be extracted.
   * @return the allowed values for the specified type (null if does not look sufficiently like an Enum).
   */
  static String[] values(final Class<?> type) {
    if (type == null) {
      return null;
    }
    if (!isValidEnum(type)) {
      return null;
    }
    try {
      final String valuesMethod = "values";
      final Method m = type.getMethod(valuesMethod);
      final Class<?> returnType = m.getReturnType();
      if (returnType.isArray()) {
        final Object[] ret = (Object[]) m.invoke(null);
        final String[] res = new String[ret.length];

        for (int i = 0; i < ret.length; i++) {
          res[i] = ret[i].toString().toLowerCase(Locale.getDefault()); // List enums as lowercase by default
        }
        return res;
      }
      return null;
    } catch (final NoSuchMethodException | IllegalAccessException | InvocationTargetException e) {
      // Should never happen
      throw new RuntimeException(e);
    }

  }

  /**
   * Check if looks sufficiently like an Enum to be treated as one.
   * Must implement
   *    static T[] values()
   *    static T value(String)
   * @param type class type
   * @return true iff is an Enum or looks sufficiently like one.
   */
  static boolean isValidEnum(final Class<?> type) {
    if (type == null) {
      return false;
    }
    if (type.isEnum()) {
      return true;
    }

    final Method m;
    try {
      final String valuesMethod = "values";
      m = type.getDeclaredMethod(valuesMethod);
      if (m == null) {
        return false;
      }
    } catch (final SecurityException e) {
      return false;
    } catch (final NoSuchMethodException e) {
      return false;
    }
    if (!isStatic(m)) {
      return false;
    }
    final Class<?> returnType = m.getReturnType();
    if (!returnType.isArray()) {
      return false;
    }

    final Method v;
    try {
      final String valueOfMethod = "valueOf";
      v = type.getMethod(valueOfMethod, String.class);
      if (v == null) {
        return false;
      }

    } catch (final SecurityException e) {
      return false;
    } catch (final NoSuchMethodException e) {
      return false;
    }
    if (!isStatic(v)) {
      return false;
    }
    final Class<?> returnTypev = v.getReturnType();
    if (!type.isAssignableFrom(returnTypev)) {
      return false;
    }
    return true;
  }

  static boolean isStatic(final Method method) {
    return Modifier.isStatic(method.getModifiers());
  }

  /**
   * Sets the maximum number of times the flag can be specified.
   *
   * @param count the maximum number of times the flag can be specified.
   * @return this flag, so calls can be chained.
   */
  public Flag setMaxCount(final int count) {
    if ((count < 1) || (count < mMinCount)) {
      throw new IllegalArgumentException("MaxCount (" + count
          + ") must not be 0 or less than MinCount (" + mMinCount + ")");
    }
    mMaxCount = count;
    return this;
  }

  /**
   * Gets the maximum number of times the flag can be specified.
   *
   * @return the maximum number of times the flag can be specified.
   */
  public int getMaxCount() {
    return mMaxCount;
  }

  /**
   * Sets the minimum number of times the flag can be specified.
   *
   * @param count the minimum number of times the flag can be specified.
   * @return this flag, so calls can be chained.
   */
  public Flag setMinCount(final int count) {
    if (count > mMaxCount) {
      throw new IllegalArgumentException("MinCount (" + count
          + ") must not be greater than MaxCount (" + mMaxCount + ")");
    }
    if (count == Integer.MAX_VALUE) {
      throw new IllegalArgumentException(
      "You're crazy man -- MinCount cannot be Integer.MAX_VALUE");
    }
    mMinCount = count;
    return this;
  }

  /**
   * Gets the minimum number of times the flag can be specified.
   *
   * @return the minimum number of times the flag can be specified.
   */
  public int getMinCount() {
    return mMinCount;
  }

  /**
   * Return the number of times the flag has been set.
   *
   * @return the number of times the flag has been set.
   */
  public int getCount() {
    return mParameter.size();
  }

  /**
   * Return true if the flag has been set.
   *
   * @return true if the flag has been set.
   */
  public boolean isSet() {
    return mParameter.size() > 0;
  }

  /**
   * Gets the character name of this flag, if set.
   *
   * @return the character name of this flag, or null if no character name has
   * been set.
   */
  public Character getChar() {
    return mFlagChar;
  }

  /**
   * Gets the name of the flag.
   *
   * @return the name of the flag.
   */
  public String getName() {
    return mFlagName;
  }

  /**
   * Gets the description of the flag's purpose.
   *
   * @return the description.
   */
  public String getDescription() {
    return mFlagDescription;
  }

  /**
   * Gets the description of the flag parameter. This is usually a single word
   * that indicates a little more than the parameter type.
   *
   * @return the parameter description, or null for untyped flags.
   */
  public String getParameterDescription() {
    return mParameterDescription;
  }

  /**
   * Gets the type of the parameter. This will return null for untyped
   * (switch) flags. Parameters will be checked that they are of the specified
   * type.
   *
   * @return the parameter type, or null if the flag is untyped.
   */
  public Class<?> getParameterType() {
    return mParameterType;
  }

  /**
   * Gets the default value of the parameter.
   *
   * @return the default value, or null if no default has been specified.
   */
  public Object getParameterDefault() {
    return mParameterDefault;
  }

  /**
   * Sets the default value of the parameter.
   *
   * @param paramDefault a default value that can be used for optional flags.
   * @return this flag, so calls can be chained.
   */
  public Flag setParameterDefault(final Object paramDefault) {
    if (mParameterType == null) {
      throw new IllegalArgumentException("Cannot set default parameter for untyped flags");
    }
    mParameterDefault = paramDefault;
    return this;
  }

  /**
   * Defines the set of strings that are valid for this flag.
   *
   * @param range a collection of Strings.
   * @return this flag, so calls can be chained.
   */
  public Flag setParameterRange(final Collection<String> range) {
    //System.err.println("setParameterRange range=" + range.toString());
    final String[] rarray = range.toArray(new String[range.size()]);
    return setParameterRange(rarray);
  }

  /**
   * Defines the set of strings that are valid for this flag.
   *
   * @param range an array of Strings.
   * @return this flag, so calls can be chained.
   */
  public Flag setParameterRange(final String[] range) {
    if (mParameterType == null) {
      throw new IllegalArgumentException("Cannot set parameter range for no-arg flags.");
    }
    if (range == null) {
      mParameterRange = null;
    } else {
      if (range.length < 1) {
        throw new IllegalArgumentException("Must specify at least one value in parameter range.");
      }
      final List<String> l = new ArrayList<>();
      for (final String s : range) {
        try {
          parseValue(s);
        } catch (final Exception e) {
          throw new IllegalArgumentException("Range value " + s + " could not be parsed.");
        }
        l.add(s);
      }
      mParameterRange = Collections.unmodifiableList(l);
    }
    return this;
  }

  /**
   * Override the default minimum - maximum string with one representing the given range.
   * @param min the minimum
   * @param max the maximum
   */
  public void setPsuedoMinMaxRangeString(final int min, final int max) {
    final String str = minMaxUsage(min, max);
    if (str.length() == 0) {
      mPsuedoMinMaxString = null;
    } else {
      mPsuedoMinMaxString = str;
    }
  }

  /**
   * Gets the list of valid parameter values, if these have been specified.
   *
   * @return a <code>List</code> containing the permitted values, or null if
   * this has not been set.
   */
  public List<String> getParameterRange() {
    return mParameterRange;
  }

  /**
   * Get the value for this flag. If the flag was not user-set, then the
   * default value is returned (if defined). The value will have been checked
   * to comply with any parameter typing. If called on an untyped flag, this
   * will return Boolean.TRUE or Boolean.FALSE appropriately.
   *
   * @return a value for this flag.
   */
  public Object getValue() {
    return (isSet()) ? mParameter.get(0) : (mParameterType == null) ? Boolean.FALSE
        : mParameterDefault;
  }

  /**
   * Get a collection of all values set for this flag. This is for flags that
   * can be set multiple times. If the flag was not user-set, then the
   * collection contains only the default value (if defined).
   *
   * @return a <code>Collection</code> of the supplied values.
   */
  public List<Object> getValues() {
    final List<Object> result;
    if (isSet()) {
      result = mParameter;
    } else {
      result = new ArrayList<>();
      if (mParameterType == null) {
        result.add(Boolean.FALSE);
      } else if (mParameterDefault != null) {
        result.add(mParameterDefault);
      }
    }
    //System.err.println(mFlagName + ":" + result);
    return Collections.unmodifiableList(result);
  }

  void reset() {
    //System.err.println("reset");
    mParameter = new ArrayList<>();
  }

  FlagValue setValue(final String valueStr) {
    if (mParameter.size() >= mMaxCount) {
      throw new FlagCountException("Value cannot be set more than " + mMaxCount + "times for flag: " + mFlagName);
    }
    if (mParameterRange != null) {
      if (mRangeList) {
        final String[] vs = StringUtils.split(valueStr, ',');
        for (String vi : vs) {
          if (!mParameterRange.contains(vi)) {
            throw new IllegalArgumentException("A value supplied is not in the set of allowed values.");
          }
        }
      } else {
        if (!mParameterRange.contains(valueStr)) {
          throw new IllegalArgumentException("Value supplied is not in the set of allowed values.");
        }
      }
    }
    if (mRangeList) {
      final List<Object> values = new ArrayList<>();
      final String[] valueStrs = StringUtils.split(valueStr, ',');
      for (final String valueStr2 : valueStrs) {
        final Object value = parseValue(valueStr2);
        mParameter.add(value);
        values.add(value);
      }
      return new FlagValue(this, values);
    } else {
      final Object value = parseValue(valueStr);
      mParameter.add(value);
      return new FlagValue(this, value);
    }
  }

  /**
   * Converts the string representation of a parameter value into the
   * appropriate Object. This default implementation knows how to convert
   * based on the parameter type for several common types. Override for custom
   * parsing.
   */
  Object parseValue(final String valueStr) {
    return mParameterType == null ? Boolean.TRUE : Flag.instanceHelper(mParameterType, valueStr);
  }

  @Override
  public boolean equals(final Object other) {
    return other instanceof Flag && getName().equals(((Flag) other).getName());
  }

  @Override
  public int hashCode() {
    return getName() == null ? 0 : getName().hashCode();
  }

  @Override
  public int compareTo(final Flag other) {
    if (other == null) {
      return -1;
    }
    if (other.getName() != null) {
      return getName().compareTo(other.getName());
    }
    return -1;
  }

  /** Make a compact usage string (prefers char name if present). */
  String getCompactFlagUsage() {
    final StringBuilder sb = new StringBuilder();
    if (getChar() != null) {
      sb.append(CFlags.SHORT_FLAG_PREFIX).append(getChar());
    } else {
      sb.append(CFlags.LONG_FLAG_PREFIX).append(getName());
    }
    final String usage = getParameterDescription();
    if (usage.length() > 0) {
      sb.append(' ').append(usage);
    }
    return sb.toString();
  }

  /** Make a usage string. */
  String getFlagUsage() {
    final StringBuilder sb = new StringBuilder();
    sb.append(CFlags.LONG_FLAG_PREFIX).append(getName());
    if (getParameterType() != null) {
      sb.append('=').append(getParameterDescription());
    }
    return sb.toString();
  }

  static String minMaxUsage(int min, int max) {
    final StringBuilder ret = new StringBuilder();
    if (min >= 1 && max > 1) {
      if (max == Integer.MAX_VALUE) {
        ret.append("Must be specified ").append(min).append(" or more times");
      } else if (max - min == 0) {
        ret.append("Must be specified ").append(min).append(" times");
      } else if (max - min == 1) {
        ret.append("Must be specified ").append(min).append(" or ").append(max).append(" times");
      } else {
        ret.append("Must be specified ").append(min).append(" to ").append(max).append(" times");
      }
    } else {
      if (min == 0) {
        if (max > 1) {
          if (max == Integer.MAX_VALUE) {
            ret.append("May be specified 0 or more times");
          } else {
            ret.append("May be specified up to ").append(max).append(" times");
          }
        }
      }
    }
    return ret.toString();
  }

  void appendLongFlagUsage(final WrappingStringBuilder wb, final int longestUsageLength, final boolean showExtended) {
    if (!CFlags.displayFlag(showExtended, this)) {
      return;
    }

    wb.append("  ");
    if (getChar() == null) {
      wb.append("    ");
    } else {
      wb.append(CFlags.SHORT_FLAG_PREFIX).append(getChar()).append(", ");
    }

    final String usageStr = getFlagUsage();
    wb.append(getFlagUsage());
    for (int i = 0; i < longestUsageLength - usageStr.length(); i++) {
      wb.append(" ");
    }
    wb.append(" ");

    final StringBuilder description = new StringBuilder(getDescription());

    final List<String> range = getParameterRange();
    if (range != null) {
      if (mRangeList) {
        description.append(" (Must be one or more of ").append(Arrays.toString(range.toArray())).append(" in a comma separated list)");
      } else {
        description.append(" (Must be one of ").append(Arrays.toString(range.toArray())).append(")");
      }
    }
    final String minMaxUsage;
    if (mPsuedoMinMaxString != null) {
      minMaxUsage = mPsuedoMinMaxString;
    } else {
      minMaxUsage = minMaxUsage(getMinCount(), getMaxCount());
    }
    if (minMaxUsage.length() > 0) {
      description.append(". ").append(minMaxUsage);
    }

    final Object def = getParameterDefault();
    if (def != null) {
      final String defs;
      if (def instanceof Double) {
        defs = Utils.realFormat((Double) def);
      } else if (isValidEnum(mParameterType)) {
        defs = def.toString().toLowerCase(Locale.getDefault());
      } else {
        defs = def.toString();
      }
      description.append(" (Default is ").append(defs).append(")");
    }

    wb.wrapText(description.toString());
    wb.append(CFlags.LS);
  }

  private static String autoDescription(final Class<?> type) {
    final String result = type.getName();
    return result.substring(result.lastIndexOf('.') + 1).toUpperCase(Locale.getDefault());
  }

  private static final Set<String> BOOLEAN_AFFIRMATIVE = new HashSet<>();
  private static final Set<String> BOOLEAN_NEGATIVE = new HashSet<>();
  static {
    BOOLEAN_AFFIRMATIVE.add("true");
    BOOLEAN_AFFIRMATIVE.add("yes");
    BOOLEAN_AFFIRMATIVE.add("y");
    BOOLEAN_AFFIRMATIVE.add("t");
    BOOLEAN_AFFIRMATIVE.add("1");
    BOOLEAN_AFFIRMATIVE.add("on");
    BOOLEAN_AFFIRMATIVE.add("aye");
    BOOLEAN_AFFIRMATIVE.add("hai");
    BOOLEAN_AFFIRMATIVE.add("ja");
    BOOLEAN_AFFIRMATIVE.add("da");
    BOOLEAN_AFFIRMATIVE.add("ya");
    BOOLEAN_AFFIRMATIVE.add("positive");
    BOOLEAN_AFFIRMATIVE.add("fer-shure");
    BOOLEAN_AFFIRMATIVE.add("totally");
    BOOLEAN_AFFIRMATIVE.add("affirmative");
    BOOLEAN_AFFIRMATIVE.add("+5v");

    BOOLEAN_NEGATIVE.add("false");
    BOOLEAN_NEGATIVE.add("no");
    BOOLEAN_NEGATIVE.add("n");
    BOOLEAN_NEGATIVE.add("f");
    BOOLEAN_NEGATIVE.add("0");
    BOOLEAN_NEGATIVE.add("off");
  }

  static Object instanceHelper(final Class<?> type, final String stringRep) {
    try {
      if (type == Boolean.class) {
        final String lStr = stringRep.toLowerCase(Locale.getDefault());
        if (BOOLEAN_AFFIRMATIVE.contains(lStr)) {
          return Boolean.TRUE;
        } else if (BOOLEAN_NEGATIVE.contains(lStr)) {
          return Boolean.FALSE;
        } else {
          throw new IllegalArgumentException("Invalid boolean value " + stringRep);
        }
      } else if (type == Byte.class) {
        return Byte.valueOf(stringRep);
      } else if (type == Character.class) {
        return stringRep.charAt(0);
      } else if (type == Float.class) {
        return Float.valueOf(stringRep);
      } else if (type == Double.class) {
        return Double.valueOf(stringRep);
      } else if (type == Integer.class) {
        return Integer.valueOf(stringRep);
      } else if (type == Long.class) {
        return Long.valueOf(stringRep);
      } else if (type == Short.class) {
        return Short.valueOf(stringRep);
      } else if (type == File.class) {
        return new File(stringRep);
      } else if (type == URL.class) {
        return new URL(stringRep);
      } else if (type == String.class) {
        return stringRep;
      } else if (isValidEnum(type)) {
        return valueOf(type, stringRep.toUpperCase(Locale.getDefault()));
      } else if (type == Class.class) {
        return Class.forName(stringRep);
      } else if (type == IntegerOrPercentage.class) {
        return IntegerOrPercentage.valueOf(stringRep);
      }
    } catch (final MalformedURLException e) {
      throw new IllegalArgumentException("Badly formatted URL: " + stringRep);
    } catch (final NumberFormatException e) {
      throw new IllegalArgumentException("");
    } catch (final ClassNotFoundException e) {
      throw new IllegalArgumentException("Class not found: " + stringRep);
    }
    throw new IllegalArgumentException("Unknown parameter type: " + type);
  }

  /**
   * When set to true, this flag can take a comma-separated list of range values
   * and produce an list of those values
   * @param rangeList true if flag should represent a list of the range values, false otherwise.
   * @return this flag, so calls can be chained.
   */
  public Flag setRangeList(final boolean rangeList) {
    mRangeList = rangeList;
    return this;
  }

  /**
   * @param category the category to set
   * @return this flag, so calls can be chained.
   */
  public Flag setCategory(final String category) {
    mCategory = category;
    return this;
  }

  /**
   * @return the category
   */
  public String getCategory() {
    return mCategory;
  }
}
