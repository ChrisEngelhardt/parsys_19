// https://github.com/chapel-lang/chapel/blob/master/test/exercises/MonteCarloPi/deitz/MonteCarloPi/multiLocaleTaskParallelMonteCarloPi.chpl
use Random;

config const n = 1000000000;
config const tasks = 4; //replace with number of tasks needed
config const seed = 1231238923;

writeln("Number of locales   = ", numLocales);
writeln("Number of points    = ", n);
writeln("Random number seed  = ", seed);
writeln("Number of tasks     = ", tasks, " (per locale)");

//
// On each locale, for each task on that locale, construct a
// RandomStream object, run the Monte Carlo simulation, and delete the
// RandomStream object.  Store the resulting count in an array of
// counts, one element per task per locale.  Since there are no
// parallel accesses to the RandomStream object (each task has its own
// object), set parSafe to false to avoid locking overhead.
//
var counts: [LocaleSpace] [1..tasks] int;
coforall loc in Locales do on loc {
  var myN = (loc.id+1)*n/numLocales - (loc.id)*n/numLocales;
  coforall task in 1..tasks {
    var rs = new borrowed NPBRandomStream(real, seed + loc.id*tasks*2 + task*2, parSafe=false);
    var count = 0;
    for i in (task-1)*myN/tasks+1..task*myN/tasks do
      count += rs.getNext()**2 + rs.getNext()**2 <= 1.0;
    counts[loc.id][task] = count;
  }
}

var count = 0;
for loc in Locales do
  for task in 1..tasks do
    count += counts[loc.id][task];

writef("Approximation of PI = %{#.#######}\n", count * 4.0 / n);
