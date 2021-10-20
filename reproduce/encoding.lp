% ------------------------------------------------------------------------------
% This encoding tracks the arrival at a vertex but also the exit time from that
% vertex.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Assign tasks to robots.
% ------------------------------------------------------------------------------
% #include "common/task_assignment.lp".

1 { assignment(R,T) : robot(R) } 1 :- task(T,_), not depends(deliver,_,T).
assignment(R,T')                   :- assignment(R,T), depends(deliver,T,T').

same_robot(T,T')  :- assignment(R,T), assignment(R,T'), T < T',
                     not depends(deliver,T,T').
same_robot(T,T')  :- depends(deliver,T,T').


% ------------------------------------------------------------------------------
% Task sequencing - don't choose pickup-putdown sequences since they are
% completely determined by the pickup.
% ------------------------------------------------------------------------------
% #include "common/task_sequencing_direct.lp".

0 { task_sequence(T,T') : depends(deliver,T',_) } 1 :- depends(deliver,_,T).
task_sequence(T,T')                                 :- depends(deliver,T,T').

:- task_sequence(T,T'),
   not same_robot(T,T'), not same_robot(T',T), not depends(deliver,T,T').

% Make sure a task can't have more than one predecessor.

:- task(T,_), 2 #count{ T' : task_sequence(T',T) }.

% Each assigned robot must have exactly one linear sequence; hence must have
% only a single first task.  Note: could have some cyclic sequences as well but
% we don't prevent that here.

:- assignment(R,_), not #count{ T : assignment(R,T), not task_sequence(_,T) } = 1.


% ------------------------------------------------------------------------------
% Map each task, task assignment, and task sequence to paths.
% ------------------------------------------------------------------------------
% #include "path/task_to_path.lp".

path(T,V)             :- task(T,V).
path(R,V)             :- assignment(R,_), home(R,V).
path(R,V)             :- home(R,V), not start(R,V).

path_sequence(T,T')   :- task_sequence(T,T').
path_sequence(T,R)    :- task_sequence(_,T), not task_sequence(T,_), assignment(R,T).

path_assignment(R,T)  :- assignment(R,T).
path_assignment(R,R)  :- assignment(R,_).
path_assignment(R,R)  :- home(R,V), not start(R,V).


% ------------------------------------------------------------------------------
% Generate a set of moves for each path. The path is acyclic and the start and
% end must match the expected start and end positions to satisfy the task.
% ------------------------------------------------------------------------------
% #include "path/move_choice_domain.lp".

0 { move(T,V,V') : edge(V,V',_) } 1 :- task(T,_), edge(V,_,_).
0 { move(T,V,V') : edge(V,V',_) } 1 :- task(T,_), edge(_,V',_).

0 { move(R,V,V') : edge(V,V',_) } 1 :- robot(R), edge(V,_,_).
0 { move(R,V,V') : edge(V,V',_) } 1 :- robot(R), edge(_,V',_).
:- robot(R), not path(R,_), move(R,_,_).


% ------------------------------------------------------------------------------
% Ensure that the move choices for a path form a proper sequence and ensure that
% each path starts and ends and the appropriate vertices.
% ------------------------------------------------------------------------------
% #include "path/move_sequencing.lp".

% ------------------------------------------------------------------------------
% Define the first and last visit based on the moves
% ------------------------------------------------------------------------------

first_visit(P,V) :- move(P,V,_), not move(P,_,V).
last_visit(P,V) :-  move(P,_,V), not move(P,V,_).

% Special case if a robot is already at its path destination
last_visit(P,V) :- path(P,V), not move(P,_,_).
first_visit(P,V) :- path(P,V), not move(P,_,_).

% ------------------------------------------------------------------------------
% Because the move choices ensure that for any vertex connected to a path there
% can be at most one incoming edge and at most one outgoing edge, we therefore
% only need to look at the last vertex visited to ensure that a path is a single
% unbroken sequence that ends at the appropriate vertex.
% ------------------------------------------------------------------------------

:- #count{ V : last_visit(P,V) } > 1, path(P,_).
:- path(P,V), not last_visit(P,V).

% ------------------------------------------------------------------------------
% Match the start and end of paths associated with the same robot
% ------------------------------------------------------------------------------

% First visit for the first path - two cases if robot has assigments or not.
:- start(R,V), path_assignment(R,P),
   path_sequence(P,_), not path_sequence(_,P), not first_visit(P,V).
:- start(R,V), path(R,_), not assignment(R,_), not first_visit(R,V).

% Path sequences
:- path_sequence(P,P'), path(P,V), not first_visit(P',V).

% ------------------------------------------------------------------------------
% For a robot with nothing to do (so it starts at its home vertex and has no
% assignments) - there should be no moves to or from its home vertex.
% ------------------------------------------------------------------------------

%:- start(R,V), not path(R,_), visit(_,V).

:- start(R,V), not path(R,_), move(_,V,_).
:- start(R,V), not path(R,_), move(_,_,V).

% ------------------------------------------------------------------------------
% Track the moves that are visited.
% ------------------------------------------------------------------------------

visit(P,V) :- move(P,V,_).
visit(P,V) :- path(P,V).

% ------------------------------------------------------------------------------
% Decide the order in which robots visit a vertex (as part of a given path) to
% ensure a conflict free ordering.
% ------------------------------------------------------------------------------
%#include "path/visit_ordering.lp".

same_walk(T,T')                :- same_robot(T,T').
same_walk(T',T)                :- same_robot(T,T').
same_walk(R,T)                 :- assignment(R,T).

% P visits V before P' visit V'
{ before((P,V),(P',V')) }    :- visit(P,V), visit(P',V'),
                                conflict(v,V,V'), P < P',
                                not same_walk(P,P').
before((P',V'),(P,V))        :- visit(P,V), visit(P',V'),
                                conflict(v,V,V'), P < P',
                                not same_walk(P,P'),
                                not before((P,V),(P',V')).

% The start of a robot's walk must be first at the starting vertex (there are
% two cases: robot with/without an assignment)
:- start(R,V), path(R,_), not assignment(R,_), before((_,_),(R,V)).
:- start(R,V), path_assignment(R,P),
   path_sequence(P,_), not path_sequence(_,P), before((_,_),(P,V)).

% The end of a robot's walk must arrive at its destination last
:- robot(R), path(R,V), before((R,V),(_,_)).

% Ensure the transition between two paths is reflected in the before/3 relation
:- before((P,V),(P',V')), path(P,V), path_sequence(P,P''),
   not before((P'',V),(P',V')).

% ------------------------------------------------------------------------------
% Track successive moves to make sure there are no swap or overtaking actions
% ------------------------------------------------------------------------------

% Prevent robots that are moving in the same direction from overtaking. If two
% robot move in the same directions over two nodes then one robot must complete
% both moves before the other.
:- move(P,V1,V2), move(P',V1',V2'),
   conflict(v,V1,V1'), conflict(v,V2,V2'), before((P,V1),(P',V1')),
   not before((P,V2),(P',V2')).


% Prevent robots that are moving in the opposite direction from swapping. If two
% robot move in opposite directions over two nodes then one robot must complete
% both moves before the other.
:- move(P,V1,V2), move(P',V1',V2'),
   conflict(v,V1,V2'), conflict(v,V2,V1'), before((P,V1),(P',V2')),
   not before((P,V2),(P',V1')).

% ------------------------------------------------------------------------------
% Travel times must satisfy edge weights. Transition between tasks must satisfy
% epsilon time delay.
% ------------------------------------------------------------------------------

#const epsilon=10.

% Travel times along a vertex and edge - note: a robot "exits" its final
% destination at the upper bound.

&diff{ arrive(P,V) - exit(P,V) } <= 0        :- visit(P,V).
&diff{ exit(P,V) - arrive(P,V') } <= -W      :- move(P,V,V'), edge(V,V',W).

% Action time
&diff{ arrive(T,V) - exit(T,V) } <= -epsilon :- visit(T,V), task(T,V).

% Match end of one task and start of next
&diff{ exit(P,V) - arrive(P',V) } <= 0       :- path_sequence(P,P'), path(P,V).


% Each robot must "arrive" at its starting location at exactly time 0. Note:
% must be exact to avoid boundary case where a robot "appears" only after time 0.
&diff{ 0 - arrive(P,V) } <= 0    :- path(P,_), not path_sequence(_,P),
                                    path_assignment(R,P), start(R,V).
&diff{ arrive(P,V) - 0 } <= 0    :- path(P,_), not path_sequence(_,P),
                                    path_assignment(R,P), start(R,V).

% Returning home before the upper bound.
&diff{ arrive(R,V) - bound } <= 0      :- home(R,V).
&diff{ exit(R,V) - bound } <= 0        :- home(R,V).
&diff{ bound - exit(R,V) } <= 0        :- home(R,V).

% ------------------------------------------------------------------------------
% Handle the before((P,V),(P',V')) predicate.  NOTE: We are modelling a vertex
% as taking up space enough to contain a robot. Therefore we can have a robot
% exit a node at the same instance that another robot arrives at that node. But
% we are also allowing a robot to arrive and exit a node at the same
% instant. This creates a problem that before((P,V),(P',V)) allows P and P' to
% both arrive and exit the node at the same instant without violating the
% constraint. We could add a delay between arriving and exiting, but instead
% make sure two robots cannot arrive at the same time.
% ------------------------------------------------------------------------------

&diff{ arrive(P,V) - arrive(P',V') } <= -1  :- before((P,V),(P',V')).
&diff{ exit(P,V) - arrive(P',V') } <= 0     :- before((P,V),(P',V')).

% ------------------------------------------------------------------------------
% Two cases for edge conflicts that are not subsumed by the vertex conflicts.
% ------------------------------------------------------------------------------

&diff{ arrive(P,V2) - exit(P',V1') } <= 0  :-
                                       before((P,V2),(P',V2')),
                                       move(P,V1,V2), move(P',V1',V2'),
                                       conflict(e,(V1,V2),(V1',V2')).

&diff{ arrive(P,V2) - exit(P',V1') } <= 0  :-
                                       before((P,V1),(P',V2')),
                                       move(P,V1,V2), move(P',V1',V2'),
                                       conflict(e,(V1,V2),(V1',V2')).

% ------------------------------------------------------------------------------
% Assignment dependency checking only required for wait dependency because
% deliver dependency is for a single robot and handled by task sequencing.
% Note: the following relies on replacement/wait dependency being at the same
% vertex so the before/3 constraints will ensure that the two robots cannot be
% on the vertex at the same time.
% ------------------------------------------------------------------------------

&diff{ exit(T,V) - arrive(T',V')} <= 0  :- depends(wait,T,T'),
                                           task(T,V), task(T',V').

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
