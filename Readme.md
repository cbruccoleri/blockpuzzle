# Tile Puzzle (15)

A tile puzzle game in which the user must sort numerical tiles in a given order.

## Description

The game is playable, but the most interesting part is probably the
implementation of search algorithms. This program is a demonstration of some
basic search methods used in classical AI. When the solution space is
enumerable, finite, and the combinatorial expansion can be kept under control,
search algorithms are a very effective method to find solutions to hard
problems. Puzzle games such as this are a great example for some of these
techniques.

This very simple game is also a demonstration of the OneLoneCoder (javidx9,
David Barr) Pixel Game Engine, a simple and versatile engine to build algorithm
visualizations, UI, and games.

## Playing

You can move the tiles with the arrow, shuffle them with the spacebar and solve
automatically with F2.

    F1: Output this help to stdout.
    F2: Execute Depth-First Search on the current configuration.
    F3: Replay found solution.
    F4: Swap two random tiles.
Arrows: Move the blank tile manually
 Space: Shuffle the tiles

## Implementation

Since I decided to use C++ for performance, I also opted for using modern C++
and made liberal use of the Standard Template Library (STL), smart pointers,
etc. The same algorithms implemented in Python would look quite a bit simpler,
but would run slower. When one is writing algorithms on NP-hard problems, speed
matters.

There are many improvements to be made, both in style and efficiency. Since this
is only a fun demo, I think that readability matters, so I kept it fairly
simple.

Copyright (c) 2020 Christian Bruccoleri - LICENSE: MIT

The Pixel Game Engine is licensed under its own terms (see OLC-3 in the header
file).

Version: 0.1
Date:    2020-08-01

## TODO

Implement solution methods with other classical search algorithms, such as
breadth-search, Dijkistra, and A*.

- Implementing A* is the highest priority.