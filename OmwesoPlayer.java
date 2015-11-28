package s260534500;

/*RULES:
 * 32 SEEDS, 2X8=16 PITS
 * LEGALMOVE= MORE THAN 1 SEED
 * COUNTER CLOCKWISE MOVE, 1 SEED IN EACH PIT
 * END IN:
 * 1-EMPTY: TURN OVER
 * 2-OCCUPIED & 2ND ROW: OPPONENT BOTH OCCUPIED, CAPTURE; 
 * 		SOWING FROM WHERE 1ST SOWING BEGAN
 * 3-OCCUPIED BUT ELSE: TAKE ALL SEEDS, SOWING FROM THERE 
 * END/WIN:
 * OPPONENT HAS ONLY 0 OR 1 SEED IN ALL HIS PITS
 *
 *NOTE TO SELF
 * 
 * THINGS TO THINK ABOUT:
 * OPPONENT REACTION
 * STARTING STRATEGY
 * FUNCTIONS:
 * 		CHECK FOR DOUBLE OCCUPIED
 * A HEURISTIC FOR THE PRUNING
 * HASH TABLE?!?!
 * WHICH ALGO TO USE
 * 		ALPHA-BETA PRUNING WITH POSITION COST AND ACTION COST
 * OPENING
 * 		ALPHA-BETA PRUNING TO THE 5ND LEVEL
 * MIDDLEGAME(REACH EITHER 16 OR 48 SEEDS...#TBD)
 * 		MONTE CARLO
 * ENDGAME(LOSING OR WINING)
 * 		DON'T THINK I WILL CARE ABOUT THIS PART
 * 
 * REALLY GREEDY WAY:
 * -CAPTURE/SOWING THE MOST SEEDS
 * -CAPTURE
 * -TRY TO RELEASE UR MOST OCCUPIED PIT
 * -THE LEAST DAMAGING MOVE
 * 			CREATE THE LEAST VULNERABLE PIT POSITION(NO DOUBLE OCCUPIED)
 * 			IF EVER VERLNERBLE, NOT ON THE MOST OCCUPIED
 */

import boardgame.Board;
import boardgame.BoardState;
import boardgame.Move;
import boardgame.Player;
import omweso.CCBoardState;
import omweso.CCBoard;
import omweso.CCMove;

import java.util.ArrayList;
import java.util.Random;

import s260534500.mytools.MyTools;
import s260534500.mytools.Node;

public class s260534500Player extends Player {
    Random rand = new Random();

    /** You must provide a default public constructor like this,
     * which does nothing but call the base-class constructor with
     * your student number. */
    public s260534500Player() { super("260534500"); }
    public s260534500Player(String s) { super(s); }

    /** Leave this method unchanged. */
    public Board createBoard() { return new CCBoard(); }

    /** Use this method to take actions when the game is over. */
    public void gameOver( String msg, BoardState bs) {
        CCBoardState board_state = (CCBoardState) bs;

        if(board_state.haveWon()){
            System.out.println("I won!");
        }else if(board_state.haveLost()){
            System.out.println("I lost!");
        }else if(board_state.tieGame()){
            System.out.println("Draw!");
        }else{
            System.out.println("Undecided!");
        }
    }

    public Move chooseMove(BoardState bs)
    {
        // Cast the arguments to the objects we want to work with.
        CCBoardState board_state = (CCBoardState) bs;
        
        double alpha=Double.NEGATIVE_INFINITY;
        double beta=Double.POSITIVE_INFINITY;

        if(!board_state.isInitialized()){
            // Code for setting up our initial position. Also, if your agent needs to be
            // set-up in any way (e.g. loading data from files) you should do it here.

            //CCBoardState.SIZE is the width of the board, in pits.
            int[] initial_pits = new int[2 * CCBoardState.SIZE];

            // Make sure your initialization move includes the proper number of seeds.
            int num_seeds = CCBoardState.NUM_INITIAL_SEEDS;

            if(board_state.playFirst()){
                // If we're going to play first in this game
            	MyTools.setup2(initial_pits);	
            }else{
                // If we're going to play second this game
            	MyTools.setup2(initial_pits);
            }

            return new CCMove(initial_pits);
        }else{
            // Alpha-beta pruning
            Node bestmove;
            //ArrayList<CCMove> moves = board_state.getLegalMoves();
            
            //check till depth 3
            bestmove=MyTools.alphabeta(board_state, 5, alpha, beta, 0);
            return bestmove.getMove();
        }
    }
}
