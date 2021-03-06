#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


struct body{
	double x;
	double y;
	double vx;
	double vy;
	double fx;
	double fy;
	double m;
};


#define dT 1.0f		//timestep
#define G 1.0f		//Gravity constant
#define FS 100		//Field size
#define IS 1.0f		//Initial max speed

int N = 100;		//Number of bodies
int S = 10000;		//Number of steps

float randF( float MaxVal ){
	return ( (float)rand( ) / (float)RAND_MAX ) * MaxVal;
}


void sampleBodies(struct body* bodies){
	for (int i = 0;i<N;i++) {
		bodies[i].x = randF(FS);
		bodies[i].y = randF(FS);
		bodies[i].vx = randF(IS)*2-IS;
		bodies[i].vy = randF(IS)*2-IS;
		bodies[i].fx = 0;
		bodies[i].fy = 0;
		bodies[i].m = randF(1);
	}
}

void simulate(struct body* bodies){
	//calculate forces
	for (int i=0;i<N;i++) {
		struct body* b1 = &bodies[i];
		b1->fx = 0.0;
		b1->fy = 0.0;
		for (int j=0;j<N;j++) {
			if (i == j) continue;
			struct body* b2 = &bodies[j];
			double dx = b2->x-b1->x;
    		double dy = b2->y-b1->y;
			double r = sqrt(dx*dx+dy*dy);
			double force = G*(b1->m*b2->m)/(r*r);
			b1->fx += force * dx / r;
			b1->fy += force * dy / r;
		}
	}

	//update velocity/position
	for (int i=0;i<N;i++){
		struct body* b1 = &bodies[i];
		b1->vx += dT*b1->fx/b1->m;
		b1->vy += dT*b1->fy/b1->m;
		b1->x += b1->vx*dT;
		b1->y += b1->vy*dT;
	}
}

int main(int argc, char *argv[]) {
	clock_t begin = clock();

	srand((unsigned int)time(NULL));
	if (argc > 1) N = atoi(argv[1]);
	if (argc > 2) S = atoi(argv[2]);

	struct body * bodies = malloc(sizeof(struct body)*N);

	sampleBodies(bodies);

	for (int i = 0;i<S;i++) {
		simulate(bodies);
	}
	free(bodies);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%f\n",time_spent);
}