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
package com.rtg.util.memory;

import java.io.Serializable;
import java.lang.ref.Reference;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;

import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;

/**
 * Computes statistics about memory which can be seen from a particular
 * object. <p>
 *
 * Can be used two different ways: <p>
 *
 * The simple way:<br>
 * <code>MemoryUsage.size(obj)</code>
 * returns the total size (including soft references)
 * accessible from <code>obj</code>. <p>
 *
 * The complicated way: <code>MemoryUsage mem = new MemoryUsage(obj);</code>
 * <code>System.out.println(mem);</code> This prints a nice formatted
 * report of memory usage including a breakdown by class and hard
 * versus soft references.
 *
 */
public class MemoryUsage {

  /** Boundary to which objects are rounded. */
  public static final int ROUND_SIZE = 8;
  /** Minimum size of an object. */
  public static final int OBJECT_SIZE;
  /** Minimum size of an array. */
  public static final int ARRAY_SIZE = 12; //but needs to be rounded up
  /** Size of a reference. */
  public static final int REF_SIZE = 8; //4 for 32-bit JVMs
  /** Additional size of an inner class. */
  public static final int INNER_SIZE = 4;
  /** The size of one boolean as a field. */
  public static final int BOOLEAN_SIZE = 4;
  /** The size of one byte as a field. */
  public static final int BYTE_SIZE = 4;
  /** The size of one char as a field. */
  public static final int CHAR_SIZE = 4;
  /** The size of one short as a field. */
  public static final int SHORT_SIZE = 4;
  /** The size of one int as a field. */
  public static final int INT_SIZE = 4;
  /** The size of one float as a field. */
  public static final int FLOAT_SIZE = 4;
  /** The size of one long as a field. */
  public static final int LONG_SIZE = 8;
  /** The size of one double as a field. */
  public static final int DOUBLE_SIZE = 8;
  /** The size of one boolean when inside an array. */
  public static final int BOOLEAN_ARRAY = 1;
  /** The size of one byte when inside an array. */
  public static final int BYTE_ARRAY = 1;
  /** The size of one char when inside an array. */
  public static final int CHAR_ARRAY = 2;
  /** The size of one short when inside an array. */
  public static final int SHORT_ARRAY = 2;
  /** The size of one int when inside an array. */
  public static final int INT_ARRAY = 4;
  /** The size of one float when inside an array. */
  public static final int FLOAT_ARRAY = 4;
  /** The size of one long when inside an array. */
  public static final int LONG_ARRAY = 8;
  /** The size of one double when inside an array. */
  public static final int DOUBLE_ARRAY = 8;
  /** All objects encountered via hard references. Set&lt;Object&gt; */
  protected transient IdentitySet mHardrefs = null;
  /**
   * Objects encountered as Reference(s). They are held as Reference(s)
   * and so can potentially be garbage collected before being looked at
   * the second time.
   */
  private IdentitySet mSoftrefs = null;
  /** Contains objects to be excluded from scans. */
  private IdentitySet mExcludeObjects = null;
  /** Contains types to be excluded from scans. */
  private IdentitySet mExcludeClasses = null;
  /** Calculates totals for hard references. */
  TotalAndClasses mHardRefs = null;
  /** Calculates totals for soft references. */
  TotalAndClasses mSoftRefs = null;
  /** Total of memory usage across both hard and soft references. */
  private Count mTotalCount = null;
  /** The object that has been analyzed. */
  private Object mObj;

  private static final FormatReal PERCENTAGE = new FormatReal(3, 2);
  private static final FormatReal AVG = new FormatReal(8, 1);
  private static final FormatInteger BYTES = new FormatInteger(11, true);
  private static final FormatInteger COUNT = new FormatInteger(7, true);

  /**
   * Maintains memory statistics. The total number of instances of
   * objects (including arrays) and the total number of bytes consumed.
   * Used both for totals and for class specific information.
   *
   */
  public static class Count implements Serializable {
    /** Bytes of memory. */
    long mTotal = 0;
    /** Bytes of memory including all descendents. */
    long mReturn = 0;
    /** Number of objects. */
    long mCount = 0;

    /**
     * Keeps track of the depth of recursion. Used to control when adds
     * are done to <code>mReturn</code>.
     */
    private int mDepth = 0;

    /**
     * Get the total memory used in bytes.
     * @return the total memory used in bytes.
     */
    public long getBytes() {
      return mTotal;
    }

    /**
     * Get the total memory in all descendents.
     * @return the total memory used in bytes.
     */
    public long getReturned() {
      return mReturn;
    }

    /**
     * Get the total number of objects (including arrays).
     * @return the total number of objects (including arrays).
     */
    public long getCount() {
      return mCount;
    }

    /**
     * Increment the number of objects and the total memory.
     * @param bytes the amount to increment the number of bytes by.
     */
    void increment(final long bytes) {
      if (bytes < 0) {
        throw new IllegalArgumentException();
      }
      mCount++;
      mTotal += bytes;
    }

    /**
     * Increment the number of objects and the total memory and
     * returned memory.
     *
     * @param bytes the amount to increment the number of bytes by.
     * @param returned the amount to increment the returned bytes by.
     */
    void increment(final long bytes, final long returned) {
      if (bytes < 0) {
        throw new IllegalArgumentException();
      }
      mCount++;
      mTotal += bytes;
      if (mDepth <= 1) {
        mReturn += returned;
      }
    }

    /**
     * Increment the number of objects and the total memory.
     *
     * @param count Description.
     */
    void increment(final Count count) {
      mCount += count.mCount;
      mTotal += count.mTotal;
    }

    /** Increase recursion depth by one. */
    void push() {
      mDepth++;
    }

    /** Decrease recursion depth by one. */
    void pop() {
      mDepth--;
      if (mDepth < 0) {
        throw new IllegalStateException("Too many pops in memory calculation.");
      }
    }
  }

  /**
   * <code>TotalAndClasses</code> description here.
   */
  public final class TotalAndClasses implements Comparator<Class<?>> {

    /** Totals for all objects accessible from the root object. */
    final Count mTotal = new Count();

    /** Totals per class. SortedMap&lt;Class,Count&gt; */
    final Map<Class<?>, Count> mClassCounts;

    private final boolean mHard;

    /**
     * (Recursively) extract the memory accessible from an object.
     *
     * @param hard Description.
     */
    TotalAndClasses(final boolean hard) {
      mClassCounts = new TreeMap<>(this);
      mHard = hard;
    }

    int memSize(final Object obj) {
      if (obj == null) {
        return 0;
      }
      //System.err.println("memSize 1:"+System.identityHashCode(obj)+":"+obj);

      //see if have encountered already
      if (mHardrefs.contains(obj)) {
        return 0;
      }
      // See if the object is in the exclude list
      if (mExcludeObjects.contains(obj)) {
        return 0;
      }
      //remember we have encountered it
      //on first pass see if a Reference - if so remember and pick up on second pass.
      if (mHard && (obj instanceof Reference)) {
        mSoftrefs.add(obj);
        return 0;
      }
      mHardrefs.add(obj);
      //will contain softrefs on 2nd pass.

      //get the class info and calculate size
      final Class<?> classId = obj.getClass();
      // See if the class is in the exclude list
      if (mExcludeClasses.contains(classId)) {
        return 0;
      }
      final Count count = getCount(mClassCounts, classId);
      count.push();
      //System.err.println("memSize 2:"+classId+":"+System.identityHashCode(obj)+":"+obj);
      final ClassMemory.Info info = ClassMemory.getMemoryInfo(classId);
      final int bytes;
      int returnSize = 0;
      if (classId.isArray()) {
        final int length = Array.getLength(obj);
        final int arraySize = info.getArraySize();
        bytes = roundUp(ARRAY_SIZE + length * arraySize);
        if (info.isArray()) {
          //but not primitive array
          for (int i = 0; i < length; i++) {
            final Object inst = Array.get(obj, i);
            returnSize += memSize(inst);
          }
        }
      } else {
        //Ordinary object or primitive
        //scan all subfields containing Objects and arrays
        for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
          final Field field = iter.next();
          //
          //
          //
          //
          final Object value;
          try {
            value = field.get(obj);
          } catch (final IllegalAccessException e) {
            throw new IllegalStateException("Running code does not have access to inner members", e);
          }
          returnSize += memSize(value);
        }
        bytes = info.getSize();
      }
      returnSize += bytes;
      mTotal.increment(bytes);
      count.increment(bytes, returnSize);
      count.pop();
      return returnSize;
    }
    //memsize

    @Override
    public int compare(final Class<?> x, final Class<?> y) {
      if (x == null || y == null) {
        throw new IllegalArgumentException();
      }
      final String sx = x.getName();
      final String sy = y.getName();
      return sx.compareTo(sy);
    }
  }

  /** Default constructor with no analysis performed. */
  public MemoryUsage() { }

  /**
   * Create memory usage statistics for the given object.
   * @param obj Description.
   */
  public MemoryUsage(final Object obj) {
    calculateUsage(obj);
  }

  /** Clears all internal counts and data structures to release memory. */
  public void clear() {
    mObj = null;
    mHardRefs = null;
    mSoftRefs = null;
    mTotalCount = null;
    mHardrefs = null;
    mSoftrefs = null;
    mExcludeClasses = null;
    mExcludeObjects = null;
  }

  /**
   * Calculate memory usage statistics for the given object.
   * @param obj The Object to analyze.
   */
  public void calculateUsage(final Object obj) {
    calculateUsage(obj, new Object[0]);
  }

  /**
   * Calculate memory usage statistics for the given object.
   * @param obj The Object to analyze.
   * @param exclude these objects from the counts.
   */
  public void calculateUsage(final Object obj, final Object[] exclude) {
    calculateUsage(obj, exclude, new Class<?>[0]);
  }

  void calculateUsage(final Object obj, final Object[] excludeObjs, final Class<?>[] excludeClasses) {
    mObj = obj;
    mHardRefs = new TotalAndClasses(true);
    mSoftRefs = new TotalAndClasses(false);
    mTotalCount = new Count();
    mHardrefs = new IdentitySet();
    mSoftrefs = new IdentitySet();
    mExcludeObjects = new IdentitySet();
    mExcludeClasses = new IdentitySet();

    //Set up exclude items in list so they will be ignored
    mExcludeObjects.addAll(excludeObjs);
    mExcludeClasses.addAll(excludeClasses);

    //Do first pass creating hardrefs and softrefs.
    mHardRefs.memSize(obj);

    //Now rescan softrefs and add into memory any that still werent
    //found as hard refs after first being seen.
    for (final Iterator<Object> iter = mSoftrefs.getIterator(); iter.hasNext(); ) {
      final Reference<?> soft = (Reference<?>) iter.next();
      final Object o2 = soft.get();
      //System.err.println(soft+":"+o2);
      if (o2 != null && !mHardrefs.contains(soft)) {
        //
        mSoftRefs.memSize(soft);
      }
    }
    //release all that used memory.
    mHardrefs = null;
    mSoftrefs = null;
    mTotalCount.increment(mHardRefs.mTotal);
    mTotalCount.increment(mSoftRefs.mTotal);
  }

  /**
   * Gets the object that has been analyzed.
   * @return The object that has been analyzed.
   */
  public Object getObject() {
    return mObj;
  }

  /**
   * Round the size of an object up to the boundary the memory system uses.
   * @param size the size to be rounded up.
   * @return the rounded up size.
   */
  public static int roundUp(final int size) {
    final int rem = size % ROUND_SIZE;
    return rem > 0 ? size + ROUND_SIZE - rem : size;
  }

  /**
   * Simple version that just returns the total memory accessible from an object.
   * @param obj Description.
   * @return Description.
   */
  public static long size(final Object obj) {
    final MemoryUsage mem = new MemoryUsage(obj);
    return mem.getSize();
  }

  /**
   * Size of single instance of the class - for arrays and objects the
   * size of the reference to them.
   * @param classRef the Class being examined.
   * @return the size of an instance in an array in bytes.
   */
  static int refSize(final Class<?> classRef) {
    if (classRef.isArray() || !classRef.isPrimitive()) {
      return REF_SIZE;
    }
    //primitive and !array
    if (classRef == boolean.class) {
      return BOOLEAN_SIZE;
    }
    if (classRef == byte.class) {
      return BYTE_SIZE;
    }
    if (classRef == char.class) {
      return CHAR_SIZE;
    }
    if (classRef == short.class) {
      return SHORT_SIZE;
    }
    if (classRef == int.class) {
      return INT_SIZE;
    }
    if (classRef == float.class) {
      return FLOAT_SIZE;
    }
    if (classRef == long.class) {
      return LONG_SIZE;
    }
    if (classRef == double.class) {
      return DOUBLE_SIZE;
    }
    //
    //
    //
    //
    //
    throw new IllegalArgumentException();
  }

  /**
   * Size of single instance of the class within an array. For
   * non-primitives (array or Object) this is the size of the reference
   * to it.
   * @param classRef the Class being examined.
   * @return the size of an instance in an array in bytes.
   */
  static int arraySize(final Class<?> classRef) {
    if (classRef.isArray() || !classRef.isPrimitive()) {
      return REF_SIZE;
    }
    //primitive and !array
    if (classRef == boolean.class) {
      return BOOLEAN_ARRAY;
    }
    if (classRef == byte.class) {
      return BYTE_ARRAY;
    }
    if (classRef == char.class) {
      return CHAR_ARRAY;
    }
    if (classRef == short.class) {
      return SHORT_ARRAY;
    }
    if (classRef == int.class) {
      return INT_ARRAY;
    }
    if (classRef == float.class) {
      return FLOAT_ARRAY;
    }
    if (classRef == long.class) {
      return LONG_ARRAY;
    }
    if (classRef == double.class) {
      return DOUBLE_ARRAY;
    }
    throw new IllegalArgumentException();
  }

  /**
   * Total number of bytes used by the object.
   * @return The size value
   */
  public long getSize() {
    if (mObj == null) {
      throw new IllegalStateException("No object has been analysed");
    }
    return mTotalCount.getBytes();
  }

  /**
   * Total number of bytes accessible via hard references.
   * @return The hard size value
   */
  public long getHardSize() {
    if (mObj == null) {
      throw new IllegalStateException("No object has been analysed");
    }
    return mHardRefs.mTotal.getBytes();
  }

  /**
   * Total number of bytes only accessible via soft (also phantom and
   * weak) references.
   * @return The soft size value
   */
  public long getSoftSize() {
    if (mObj == null) {
      throw new IllegalStateException("No object has been analysed");
    }
    return mSoftRefs.mTotal.getBytes();
  }

  /**
   * An iterator over all the Map&lt;Class,Count&gt; of statistics vs Class
   * for hard references.
   * @return The hard iterator value
   */
  public Iterator<Map.Entry<Class<?>, Count>> getHardIterator() {
    if (mObj == null) {
      throw new IllegalStateException("No object has been analysed");
    }
    return mHardRefs.mClassCounts.entrySet().iterator();
  }

  /**
   * An iterator over all the Map&lt;Class&lt;?&gt;,Count&gt; of statistics vs Class
   * for memory only accessible via soft (also phantom and weak)
   * references.
   * @return The soft iterator value
   */
  public Iterator<Map.Entry<Class<?>, Count>> getSoftIterator() {
    if (mObj == null) {
      throw new IllegalStateException("No object has been analysed");
    }
    return mSoftRefs.mClassCounts.entrySet().iterator();
  }

  /**
   * Generate a nicely formatted report of all the class and soft/hard
   * statistics.
   * @return Description.
   */
  public String toString() {
    return toString(new StringBuilder()).toString();
  }

  /**
   * Generates a representation of this object to a StringBuilder.
   * @param sb the StringBuilder that will be appended to.
   * @return a reference to the StringBuilder.
   */
  public StringBuilder toString(final StringBuilder sb) {
    if (mObj != null) {
      sb.append("\nMemoryUsage for object: @");
      sb.append(System.identityHashCode(mObj));
      sb.append("  ");
      sb.append(mObj.getClass().getName());
      sb.append("\n");
      sb.append("Total size:");
      BYTES.format(sb, getSize()).append("B\n");

      sb.append("Hard size: ");
      BYTES.format(sb, getHardSize());
      sb.append("B  ");

      final double hardPerc = 100.0 * getHardSize() / getSize();
      PERCENTAGE.format(sb, hardPerc);
      sb.append("%\n");

      sb.append("Soft size: ");
      BYTES.format(sb, getSoftSize());
      sb.append("B  ");

      final double softPerc = 100.0 * getSoftSize() / getSize();
      PERCENTAGE.format(sb, softPerc);
      sb.append("%\n");

      sb.append("Hard class statistics\n");
      heading(sb);
      classToString(sb, getSize(), getHardIterator());

      sb.append("Soft class statistics\n");
      heading(sb);
      classToString(sb, getSize(), getSoftIterator());

      sb.append("End MemoryUsage\n");
    } else {
      sb.append("No object analysed.");
    }
    return sb;
  }

  void heading(final StringBuilder sb) {
    sb.append("   Count  "); //10
    sb.append("       Size  "); //13
    sb.append("     %  "); //8
    sb.append("   Average   "); //13
    sb.append("Cumulative"); //10
    sb.append("\n");
  }

  /**
   * Retrieve the memory statistics for a specified class. Used
   * internally by <code>memSize</code>. If this is the first time the class has
   * been accessed then a new entry is created.
   *
   * @param counts Map of class statistics entries.
   * @param id key to be used for retrieving the statistics.
   * @return the memory statistics.
   */
  private Count getCount(final Map<Class<?>, Count> counts, final Class<?> id) {
    final Count count = counts.get(id);
    if (count == null) {
      final Count novo = new Count();
      counts.put(id, novo);
      return novo;
    }
    return count;
  }

  private void classToString(final StringBuilder sb, final long totalSize, final Iterator<Map.Entry<Class<?>, Count>> iter) {
    while (iter.hasNext()) {
      final Map.Entry<Class<?>, Count> entry = iter.next();
      final Class<?> classId = entry.getKey();
      final Count count = entry.getValue();
      COUNT.format(sb, count.mCount).append("# ");
      BYTES.format(sb, count.mTotal).append("B ");
      final double perc = 100.0 * count.mTotal / totalSize;
      PERCENTAGE.format(sb, perc);
      sb.append("% ");

      AVG.format(sb, count.mTotal / (double) count.mCount);
      sb.append("B ");

      BYTES.format(sb, count.mReturn);
      sb.append("B ");

      sb.append(classId.getName());
      sb.append("\n");
    }
  }

  static {
    final Properties props = System.getProperties();
    final String vendor = props.getProperty("java.vm.vendor");
    if (vendor.equals("IBM Corporation")) {
      OBJECT_SIZE = 12;
    } else {
      OBJECT_SIZE = 8;
    }
  }

  /**
   * Quick test of print format. Also an example of how to use this stuff.
   * @param args The command line arguments.
   */
  public static void main(final String[] args) {
    final HashMap<String, String> map = new HashMap<>();
    //put a bit of stuff in to make it non trivial
    map.put("keya", "valuea");
    map.put("keyb", "valueb");
    map.put("keyc", "valuec");
    final MemoryUsage usage = new MemoryUsage(map);
    System.out.println(usage);
  }
}

