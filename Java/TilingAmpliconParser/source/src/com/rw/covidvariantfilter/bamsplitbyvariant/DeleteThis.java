/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.covidvariantfilter.bamsplitbyvariant;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 *
 * @author raine
 */
public class DeleteThis {
  
    public void test(Map<String,Integer> mapString, Map<Integer,String> mapInteger){
        Integer a = 4;
        a = 5;

        Integer i = mapString.get("aa");
        i = 4;
        mapString.put("aa", 4);
        String s = mapInteger.get(2);
        s = "bbbb";
        Map<Integer,Person> mp = new HashMap<>();
        mp.put(3, new Person());
        mp.get(3).age = 5;
        AtomicInteger at = new AtomicInteger(0);       
        at.incrementAndGet();
        
    }
    
    public void streamTest(List<String> stringList, Map<String,Integer> mapString){
        List<String> sl = stringList.stream()
                .map(x -> x.substring(0, 5))
                .filter(f -> f.startsWith("A"))
                .skip(1)
                .limit(5)
                .toList(); 
        List<Person> lst = new ArrayList<>();
        List<Integer> jjhjh = lst.stream().map(x -> x.age + x.b).toList();
        lst.stream().map(x -> x.age + x.b).collect(Collectors.toList());
        lst.stream().collect(Collectors.toMap(k -> k.age, v -> v.b));
       
       List<String> l = mapString.entrySet().stream().map(x -> x.getValue() + "-" + x.getKey()).toList();
       String s = mapString.entrySet().stream().map(x -> x.getValue() + "-" + x.getKey()).
               sorted((o1,o2)->o1.compareTo(o2)).
               collect(Collectors.joining(","));
       List<Person> tt = new ArrayList<>();
       Map<Integer,List<Person>> gg = tt.stream().sorted((o1,o2)->o1.age.compareTo(o2.age)).
               collect(Collectors.groupingBy(xvar -> xvar.age));
       
       for(Map.Entry<String,Integer> e : mapString.entrySet()){
           
       }
    }
    
    private class Person {
        Integer age;
        int b;
        String s;
    }
}
