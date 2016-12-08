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

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

/**
 * Provides a cache for information about the memory used by a class.
 * This is in the form of the total fixed memory used by each instance
 * of the class and a list of all fields declared for the class
 * (locally and in superclasses and irrespective of protections). Can
 * be run from command line:
 *
 * <pre>
 * java com.rtg.util.memory.ClassMemory CLASSNAME
 * </pre>
 *
 * It prints out the fields and size of the class.
 *
 */
public final class ClassMemory {

  /**
   * Holds information about memory used by a class.
   *
   */
  public static final class Info {
    /**
     * The total fixed memory (in bytes) used by the class. This
     * includes the space consumed by any super-classes. Does not
     * include any array dope vectors referenced by fields but does
     * include the minimal information used for an object, an an inner
     * object (where appropriate), or array information plus any
     * rounding up that occurs.
     */
    private final int mSize;

    /**
     * Spaced used by member variables of the class but not the basic
     * object size nor any rounding up.
     */
    final int mBaseSize;

    /**
     * Non-zero only for an array type when it is the size of a
     * reference to the component type. If the component type is an
     * object (including an array) then it is the size of a reference,
     * otherwise (for a primitive) it the number of bytes consumed by
     * each packed instance.
     */
    private final int mArraySize;

    /** The info for a super class (null if none). */
    private final Info mSuperInfo;

    /**
     * All fields in this class (excluding statics but including
     * private, transient etc.)
     */
    private final Field[] mFields;

    /**
     * All non-primitive fields in this class (excluding statics but
     * including private, transient etc.)
     */
    private final Field[] mNonPrim;

    /** The type of the class as given by the constants below. */
    private final int mType;

    static final int PRIM = 0;
    static final int PRIM_ARRAY = 1;
    static final int ARRAY = 2;
    static final int OBJECT = 3;


    Info(final int type, final int baseSize, final int size, final int arraySize, final Info parent, final Field[] fields, final Field[] nonPrim) {
      mBaseSize = baseSize;
      mSize = size;
      mArraySize = arraySize;
      mSuperInfo = parent;
      mFields = fields;
      mNonPrim = nonPrim;
      mType = type;
    }


    /**
     * The total fixed memory (in bytes) used by the class. This
     * includes the space consumed by any super-classes. Does not
     * include any array dope vectors referenced by fields but does
     * include the minimal information used for an object, an an inner
     * object (where appropriate) plus any rounding up that occurs.
     *
     * @return the fixed memory size used by an instance of an object.
     */
    public int getSize() {
      return mSize;
    }


    /**
     * Non-zero only for an array type when it is the size of a
     * reference to the component type. If the component type is an
     * object (including an array) then it is the size of a reference,
     * otherwise (for a primitive) it the number of bytes consumed by
     * each packed instance.
     *
     * @return the size used by a reference to the items in this array
     * type (otherwise 0).
     */
    public int getArraySize() {
      return mArraySize;
    }


    /**
     * Gets all fields included in a class. This includes private, but
     * not static fields, and all fields in super classes as well.
     *
     * @return an Iterator over all fields.
     */
    public Iterator<Field> getFieldIterator() {
      return new FieldIterator(this, true);
    }


    /**
     * Gets all non-primitive fields in a class. This includes private,
     * but not static fields, and all fields in super classes as well.
     *
     * @return an Iterator over non-primitive fields.
     */
    public Iterator<Field> getNonprimitiveIterator() {
      return new FieldIterator(this, false);
    }


    /**
     * Checks if a non-primitive array.
     *
     * @return true iff a non-primitive array.
     */
    public boolean isArray() {
      return mType == ARRAY;
    }


    /**
     * Checks if an array of primitives.
     *
     * @return true iff an array of primitives.
     */
    public boolean isPrimitiveArray() {
      return mType == PRIM_ARRAY;
    }


    /**
     * Checks if a primitive.
     *
     * @return true iff a primitive.
     */
    public boolean isPrimitive() {
      return mType == PRIM;
    }


    /**
     * Checks if an object.
     *
     * @return true iff an object.
     */
    public boolean isObject() {
      return mType == OBJECT;
    }


    /**
     * Need special iterator to walk up through superclass fields.
     *
     */
    private static final class FieldIterator implements Iterator<Field> {
      /** true iff primitive fields to be included. */
      private final boolean mPrim;

      /** index into current field array. */
      private int mIndex;

      /** The info object or one of the super-info objects. */
      private Info mInfo;


      FieldIterator(final Info info, final boolean prim) {
        mInfo = info;
        mIndex = 0;
        mPrim = prim;
        check();
      }


      @Override
      public boolean hasNext() {
        return !(mInfo == null);
      }


      @Override
      public Field next() {
        if (mInfo == null) {
          throw new NoSuchElementException();
        }
        final Field res = mPrim ? mInfo.mFields[mIndex] : mInfo.mNonPrim[mIndex];
        ++mIndex;
        check();
        return res;
      }


      /** Make sure index is valid if not keep moving up */
      private void check() {
        while (true) {
          final int len = mPrim ? mInfo.mFields.length : mInfo.mNonPrim.length;
          if (mIndex == len) {
            mIndex = 0;
            mInfo = mInfo.mSuperInfo;
            if (mInfo == null) {
              break;
            } else {
              continue;
            }
          } else {
            break;
          }
        }
      }


      @Override
      public void remove() {
        throw new UnsupportedOperationException();
      }
    }

    @Override
    public String toString() {
      final StringBuilder sb = new StringBuilder();
      toString(sb);
      return sb.toString();
    }


    void toString(final StringBuilder sb) {
      sb.append("size:").append(mSize).append(" bytes");
      sb.append(" array size:").append(mArraySize);
      sb.append(" type:");
      switch (mType) {
          case PRIM:
            sb.append("Primitive");
            break;
          case PRIM_ARRAY:
            sb.append("Primitive array");
            break;
          case ARRAY:
            sb.append("Object array");
            break;
          case OBJECT:
            sb.append("Object");
            break;
          default:
            throw new RuntimeException();
      }
      sb.append("\n");
      for (final Iterator<Field> iter = getFieldIterator(); iter.hasNext(); ) {
        final Field field = iter.next();
        if (!field.getType().isPrimitive()) {
          continue;
        }
        sb.append("\t");
        sb.append(field.toString());
        sb.append("\n");
      }
      sb.append("nonPrimitive fields:\n");
      for (final Iterator<Field> iter = getNonprimitiveIterator(); iter.hasNext(); ) {
        final Field field =  iter.next();
        sb.append("\t");
        sb.append(field.toString());
        sb.append("\n");
      }
    }
  }


  private static class Marker {
  }


  private static final Marker MARKER = new Marker();

  /**
   * Keeps cache of info for each class. <code>HashMap(Class,Info)</code>
   */
  private static final Map<Object, Object> INFO = new HashMap<>();

  /**
   * Private to prevent instantiation.
   */
  private ClassMemory() {
  }

  /**
   * Get a memory information object for a specified class. Uses cache
   * to prevent having to recompute it. An info object for the
   * requested class and all superclasses will be cached after the
   * call.
   *
   * @param classRef class reference
   * @return an information object for the class.
   */
  public static Info getMemoryInfo(final Class<?> classRef) {
    if (classRef == null) {
      return null;
    }
    final Object obj = INFO.get(classRef);
    if (obj == null) {
      INFO.put(classRef, MARKER);
      final Info res = makeInfo(classRef);
      INFO.put(classRef, res);
      return res;
    } else if (obj instanceof Marker) {
      return null;
    } else {
      return (Info) obj;
    }
  }


  /**
   * Make new <code>Info</code> object.
   *
   * @param classRef the class to extract <code>Info</code> from.
   * @return the new <code>Info</code> object.
   */
  private static Info makeInfo(final Class<?> classRef) {
    if (classRef == null) {
      return null;
    } else if (classRef.isArray()) {
      final Class<?> arrayClass = classRef.getComponentType();
      final int size = MemoryUsage.arraySize(arrayClass);

      final int type = arrayClass.isPrimitive() ? Info.PRIM_ARRAY : Info.ARRAY;
      //nameLen==2&&name.charAt(0)=='['?Info.PRIM_ARRAY:Info.ARRAY); //Hack hack

      //System.err.println("making array entry for:"+classRef+":"+nameLen+":"+type);
      getMemoryInfo(arrayClass);
      return new Info(type, MemoryUsage.ARRAY_SIZE, MemoryUsage.ARRAY_SIZE
        , size, null, new Field[0], new Field[0]);
    } else if (classRef.isPrimitive()) {
      return new Info(Info.PRIM, MemoryUsage.refSize(classRef), MemoryUsage.refSize(classRef)
        , 0, null, new Field[0], new Field[0]);
    } else { //an object
      final Class<?> superClass = classRef.getSuperclass();
      final Info superInfo = getMemoryInfo(superClass); //recursive call.
      final Field[] fields = classRef.getDeclaredFields();
      int len = 0;
      int lenPrim = 0;
      int total = 0;

      if (superInfo != null) {
        //System.err.println(classRef+":3+"+superInfo.mSize);
        total = superInfo.mBaseSize;
      }
      for (final Field field : fields) {
        final int modifier = field.getModifiers();
        if (!Modifier.isStatic(modifier)) {
          ++len;
          final int inc = MemoryUsage.refSize(field.getType());
          //System.err.println(i+":"+len+":"+total+":"+inc+":"+field);
          //System.err.println(classRef+":"+field+"+"+inc);
          total += inc;
          final Class<?> fieldType = field.getType();
          if (!fieldType.isPrimitive()) {
            ++lenPrim;
          }
          getMemoryInfo(fieldType);
        }
      }
      final Field[] newFields = new Field[len];
      final Field[] nonPrim = new Field[lenPrim];
      int j = 0;
      int k = 0;
      for (final Field field : fields) {
        final int modifier = field.getModifiers();
        if (!Modifier.isStatic(modifier)) {
          field.setAccessible(true);
          newFields[j] = field;
          ++j;
          if (!field.getType().isPrimitive()) {
            nonPrim[k] = field;
            ++k;
          }
        }
      }
      final int up = MemoryUsage.roundUp(total + MemoryUsage.OBJECT_SIZE);
      return new Info(Info.OBJECT, total, up, 0, superInfo, newFields, nonPrim);
    }
  }

  /**
   * <pre>
   * Usage: java com.rtg.util.memory.ClassMemory CLASS-NAME*
   * </pre>
   *
   * @param args The command line arguments.
   */
  public static void main(final String[] args) {
    int start = 0;
    boolean recursive = false;
    if (args[start].equals("-R")) {
      recursive = true;
      ++start;
    }

    for (int i = start; i < args.length; ++i) {
      System.out.println(args[i]);
      try {
        final Class<?> classRef = Class.forName(args[i]);
        final Info info = getMemoryInfo(classRef);
        if (recursive) {
          for (final Map.Entry<Object, Object> entry : INFO.entrySet()) {
            final Class<?> name = (Class<?>) entry.getKey();
            final Info infor = (Info) entry.getValue();
            System.out.println(name);
            System.out.println(infor);
            System.out.println();
          }
        } else {
          System.out.println(info);
        }
      } catch (final ClassNotFoundException e) {
        System.out.println("\tClass not found");
      }
    }
  }
}

