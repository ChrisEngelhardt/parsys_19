#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <mpi.h> 


#define G 6.67300e-11     // Gravity constant
#define FS 1.0e6          // Width of space
#define IS 1.0		      //Initial max speed
#define DELTAT 0.01       // Time increament
#define THETA 1.0         // Opening angle, for approximation in BH

typedef struct {
	double px;
	double py;
} Position;


typedef struct {
	double vx;
	double vy;
} Velocity;


typedef struct {
	double fx;
	double fy;
} Force;


typedef struct Cell  {
	int index;                    // Index into arrays to identify particle's 

	int numberSubcells;           // Number of subcells
	double mass;                  // Mass of particle of total mass of subtree
	double x, y;                  // Location of cell
	double cx, cy;                // Location of center of mass
	double width, height;         // Width and Height of cell
	struct Cell* subcells[4];     // Children
} Cell;

int N = 100;		//Number of bodies
int S = 10000;		//Number of steps
int rank; 			//Rank
int partPerRank;    // Number of particles per rank
int pindex;         // Index for rank


Position* position;   // Position of particles
Velocity* ivelocity;  // Initial velocity of particles
Velocity* velocity;   // Velocity of particles
double* mass;         // Mass of particles
Force* force;         // Force of particles
Cell* rootCell;       // Root of BH octtree


MPI_Datatype MPI_POSITION;
MPI_Datatype MPI_VELOCITY;


float randF( float MaxVal ){
	return ( rand( ) / (double)RAND_MAX ) * MaxVal;
}

void sampleBodies() {   
	for (int i = 0; i < N;i++) {
		position[i].px = randF(FS);
		position[i].py = randF(FS);
		ivelocity[i].vx = randF(IS)*2-IS;
		ivelocity[i].vy = randF(IS)*2-IS;
		mass[i] = randF(1);
	}
}

double computeDistance(Position a, Position b){
	double dx = a.px - b.px;
	double dy = a.py - b.py;
	
	return sqrt(dx * dx + dy * dy);
}



/*
 * Computes force without BH
 */
void computeForce(){
	for (int i = 0; i < partPerRank; i++) {

		force[i].fx = 0.0;
		force[i].fy = 0.0;

		for (int j = 0; j < N; j++){

			if (j == (i + pindex)){
				continue;   //Same body
			}

			double d = computeDistance(position[i + pindex], position[j]);

			double f = (G * (mass[i + pindex] * mass[j]) /  (d*d));

			force[i].fx += f * ((position[j].px - position[i + pindex].px) / d);
			force[i].fy += f * ((position[j].py - position[i + pindex].py) / d);
		}
	}
}

void computeVelocity(){
	for (int i = 0; i < partPerRank; i++) {
		velocity[i].vx += (force[i].fx / mass[i + pindex]) * DELTAT;
		velocity[i].vy += (force[i].fy / mass[i + pindex]) * DELTAT;
	}
}

void computePositions(){
	for (int i = 0; i < partPerRank; i++) {
		position[i + pindex].px += velocity[i].vx * DELTAT;
		position[i + pindex].py += velocity[i].vy * DELTAT;

		// If we go outside we change direction of speed  
		if ((position[i + pindex].px) >= FS || (position[i + pindex].px) <= 0){
			velocity[i].vx *= -1;
		}
		else if ((position[i + pindex].py >= FS) || (position[i + pindex].py) <= 0){
			velocity[i].vy *= -1;
		}
	}
}


Cell* createCell(double width, double height) {
	Cell* cell = malloc(sizeof(Cell));
	cell->mass = 0;
	cell->numberSubcells = 0;
	cell->index = -1;
	cell->cx = 0;
	cell->cy = 0;
	cell->width = width;
	cell->height = height;
	return cell;
}

void setLocationSubcells(Cell* cell, double width, double heigth){

	cell->subcells[0]->x = cell->x;
	cell->subcells[0]->y = cell->y;

	cell->subcells[1]->x = cell->x + width;
	cell->subcells[1]->y = cell->y;

	cell->subcells[2]->x = cell->x;
	cell->subcells[2]->y = cell->y + heigth;

	cell->subcells[3]->x = cell->x + width;
	cell->subcells[3]->y = cell->y + heigth;
}

void generateSubcells(Cell* cell) {
	
	// Calculate subcell dimensions
	double width  = cell->width / 2.0;
	double height = cell->height / 2.0;

	// Cell no longer a leaf
	cell->numberSubcells = 4;   
	
	// Create and initialize new subcells   
	for (int i = 0; i < cell->numberSubcells; i++) {
		cell->subcells[i] = createCell(width, height);
	}
	
	setLocationSubcells(cell, width, height);   
}


/*
 * Location of a given index within a cell
 */
int locateSubcell(Cell* cell, int index) {

	// Determine which subcell to add the body to
	if (position[index].px > cell->subcells[3]->x){
		if (position[index].py > cell->subcells[3]->y){
			return 3;
		}
		else{
			return 1;
		}
	}
	else{
		if (position[index].py > cell->subcells[3]->y){
			return 2;
		}
		else{
			return 0;
		}      
	}
}


/*
 * Adds a particle to a cell
 */
void addToCell(Cell* cell, int index) {

	if (cell->numberSubcells == 0 && cell->index == -1) {         
		cell->index = index;
		return;         
	}
			
	generateSubcells(cell);

	// The current cell's body must now be re-added to one of its subcells
	int sc1 = locateSubcell(cell, cell->index);
	cell->subcells[sc1]->index = cell->index;   

	// Locate subcell for new body
	int sc2 = locateSubcell(cell, index);
	
	cell->index = -1;

	if (sc1 == sc2)
		addToCell(cell->subcells[sc1], index);
	else 
		cell->subcells[sc2]->index = index;  
}


/*
 * Generates the octtree from all particles
 */
void generateOcttree() {
	
	// Initialize root of octtree
	rootCell = createCell(FS, FS);
	rootCell->x = 0;
	rootCell->y = 0;
	
	for (int i = 1; i < N; i++) {

		Cell* cell = rootCell;

		// Find which node to add the body to
		while (cell->numberSubcells != 0){
			int sc = locateSubcell(cell, i);
			cell = cell->subcells[sc];
		}      

		addToCell(cell, i);
	}
}


Cell* computeCellProperties(Cell* cell){
	if (cell->numberSubcells == 0) {
		if (cell->index != -1){
			cell->mass = mass[cell->index];
			cell->cx = position[cell->index].px;
			cell->cy = position[cell->index].py;
			return cell;
		}
	}
	else {   
		double tx = 0, ty = 0;
		for (int i = 0; i < cell->numberSubcells; i++) {
			Cell* temp = computeCellProperties(cell->subcells[i]);
			if (temp != NULL) {
				cell->mass += temp->mass;
				tx += cell->cx * temp->mass;
				ty += cell->cy * temp->mass;
			}
		}
		
		cell->cx = tx / cell->mass;
		cell->cy = ty / cell->mass;
	
		return cell;
	}
	return NULL;
}


void computeForceOfCell(Cell* cell, int index) {
	Position positionCell;
	positionCell.px = cell->cx;
	positionCell.py = cell->cy;
	double d = computeDistance(position[index], positionCell);

	double f = (G * (mass[index] * cell->mass) / (d*d));

	// Resolve forces in each direction
	force[index - pindex].fx += f * ((cell->cx - position[index].px) / d);
	force[index - pindex].fy += f * ((cell->cx- position[index].py) / d);
}

void computeForceofOcttree(Cell* cell, int index) {
	
	if (cell->numberSubcells == 0) {
		if (cell->index != -1 && cell->index != index) {
			computeForceOfCell(cell, index);
		}
	}
	else {
		Position positionCell;
		positionCell.px = cell->cx;
		positionCell.py = cell->cy;
		double d = computeDistance(position[index], positionCell);
		
		if (THETA > (cell->width / d)){ 
			// Use approximation
			computeForceOfCell(cell, index);         
		}
		else {
			for (int i = 0; i < cell->numberSubcells; i++) {
				computeForceofOcttree(cell->subcells[i], index);
			}
		}      
	}
}


/*
 * Computes the force with BH
 */
void computeForceBH(){
	for (int i = 0; i < partPerRank; i++) {

		force[i].fx = 0.0;
		force[i].fy = 0.0;

		computeForceofOcttree(rootCell, i + pindex);
	}
}

/*
 * Deletes the octtree
 */
void deleteOcttree(Cell* cell) {
	
	if (cell->numberSubcells == 0) {
		free(cell);
		return;
	}

	for (int i = 0; i < cell->numberSubcells; i++) {
		deleteOcttree(cell->subcells[i]);
	}

	free(cell);
}


void simulate(){  
	generateOcttree();
	computeCellProperties(rootCell);
	computeForceBH();
	deleteOcttree(rootCell);
		
	computeVelocity();
	computePositions();
	MPI_Allgather(position + (rank * partPerRank), partPerRank, MPI_POSITION, 
						position, partPerRank, MPI_POSITION, MPI_COMM_WORLD); 
}


int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
  	double start = MPI_Wtime();

	
	if (argc > 1) N = atoi(argv[1]);
	if (argc > 2) S = atoi(argv[2]);
		
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	// Create and commit new MPI Types
	MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_POSITION);
	MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_VELOCITY);
	MPI_Type_commit(&MPI_POSITION);
	MPI_Type_commit(&MPI_VELOCITY);

	// Number bodies per rank
	partPerRank = N / size;


	// Start index for current rank
	pindex = rank * partPerRank;


	// Allocate memory
	mass = (double *) malloc(N * sizeof(double));
	position = (Position *) malloc(N * sizeof(Position));
	ivelocity = (Velocity *) malloc(N * sizeof(Velocity));
	velocity = (Velocity *) malloc(partPerRank * sizeof(Velocity));
	force = (Force *) malloc(partPerRank * sizeof(Force));

	if (rank == 0){
		sampleBodies();
	}
	
	// Broadcast mass and position to all ranks
	MPI_Bcast(mass, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(position, N, MPI_POSITION, 0, MPI_COMM_WORLD);
	MPI_Scatter(ivelocity, partPerRank, MPI_VELOCITY, velocity, partPerRank, MPI_VELOCITY, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < S; i++) {
		simulate();
	}
	
	
	if (rank == 0){
		for (int i = 0; i < N; i++) {
			 printf("px=%f, py=%f\n", position[i].px, position[i].py);
		}
    	printf("%lf\n",MPI_Wtime() - start);
	}

	MPI_Finalize();

	return 0;                   
}
