<?php
/**
 * Implements the Kamada-Kawai algorithm for node layout.
 * Does not respect filter calls, and sometimes crashes when the view changes to it.
 *
 * @see "Tomihisa Kamada and Satoru Kawai: An algorithm for drawing general indirect graphs. Information Processing Letters 31(1):7-15, 1989" 
 * @see "Tomihisa Kamada: On visualization of abstract objects and relations. Ph.D. dissertation, Dept. of Information Science, Univ. of Tokyo, Dec. 1988."
 *
 * @author ggsonic
 */
class KKLayout{

	var $EPSILON = 0.1;

	var $currentIteration;
    var $maxIterations = 500;//2000;
	var $status = "KKLayout";

	var $L;			// the ideal length of an edge
	var $K = 1;		// arbitrary const number
	var $dm = array();     // distance matrix

	var $adjustForGravity = true;
	var $exchangeVertices = true;

	var $vertices=array();
	var $xydata=array();
	var $size = array('width'=>600,'height'=>600);

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    var $distance=array();

    /**
     * The diameter of the visible graph. In other words, the maximum over all pairs
     * of vertices of the length of the shortest path between a and bf the visible graph.
     */
	var $diameter = 5.0;//5.0;

    /**
     * A multiplicative factor which partly specifies the "preferred" length of an edge (L).
     */
    var $length_factor = 0.9;

    /**
     * A multiplicative factor which specifies the fraction of the graph's diameter to be 
     * used as the inter-vertex distance between disconnected vertices.
     */
    var $disconnected_multiplier = 0.7;//0.5;
    
    var $initialized = false;
    var $locations = array();
    var $graph  = array();
	/**
	 * Creates an instance for the specified graph and distance metric.
	 */
    function __construct($graph){
        $this->graph = $graph;
    }
    
    public function setSize($width,$height){
        $this->size = array('width'=>$width,'height'=>$height);
    }
    
	/**
	 * Returns true once the current iteration has passed the maximum count.
	 */
	function done() {
		if ($this->currentIteration > $this->maxIterations) {
			return true;
		}
		return false;
	}
	
    public function initialize() {
    	$this->currentIteration = 0;

    	if($this->graph != null && $this->size != null) {

        	$height = $this->size['height'];
    		$width = $this->size['width'];

    		$n = count($this->graph);
    		$this->dm = array();//new Array[n][n];
    		$this->vertices = $this->graph;
    		$this->xydata = array();

    		// assign IDs to all visible vertices
    		while(true) {
    			try {
    				$index = 0;
    				foreach($this->graph as $v) {
    					$xyd = $this->transform($v);
    					$this->vertices[$index] = $v;
    					$this->xydata[$index] = $xyd;
    					$index++;
    				}
    				break;
    			} catch(Exception $e) {}
    		}

//    		$this->diameter = 5.0;//DistanceStatistics.<V,E>diameter(graph, distance, true);

    		$L0 = min(array($height, $width));
    		$this->L = ($L0 / $this->diameter) * $this->length_factor;  // length_factor used to be hardcoded to 0.9
    		//L = 0.75 * sqrt(height * width / n);
    		for ($i = 0; $i < $n - 1; $i++) {
    			for ($j = $i + 1; $j < $n; $j++) {
    				$d_ij = $this->getDistance($this->vertices[$i], $this->vertices[$j]);
    				$d_ji = $this->getDistance($this->vertices[$j], $this->vertices[$i]);
    				$dist = $this->diameter * $this->disconnected_multiplier;
    				if ($d_ij != null)
    					$dist = min(array($d_ij, $dist));
    				if ($d_ji != null)
    					$dist = min(array($d_ji, $dist));
    				$this->dm[$i][$j] = $this->dm[$j][$i] = $dist;
    			}
    		}
    	}
	}
	
	/**
	 * init node's xydate
	 * @return ArrayObject
	 */
	protected function transform(){
	    $x = rand(10,$this->size['width']-1);
	    $y = rand(10,$this->size['height']-1);
	    return array('x'=>$x,'y'=>$y);
	}
    
	/**
	 * For now, it is just a very simple strategy to get weight of distance
	 * TODO: fix it!
	 * @param $from
	 * @param $to
	 * @return distance weight
	 */
	protected function getDistance($from, $to){
	    $fromid = $from['idx'];
	    $toid   = $to['idx'];
	    if (in_array($toid,(array)$from['conn'])){
	        $weight = 1;
	    }else{
	        $weight = null;
	    }
	    return $weight;
	}
	
	public function step() {
		try {
			$this->currentIteration++;
			$energy = $this->calcEnergy();
			$this->status = "Kamada-Kawai V=" . count($this->graph)
			. "(" + count($this->graph) + ")"
			. " IT: " + $this->currentIteration
			. " E=" + $energy
			;

			$n = count($this->graph);
			if ($n == 0)
				return;

			$maxDeltaM = 0;
			$pm = -1;            // the node having max deltaM
			for ($i = 0; $i < $n; $i++) {
				$deltam = $this->calcDeltaM($i);

				if ($maxDeltaM < $deltam) {
					$maxDeltaM = $deltam;
					$pm = $i;
				}
			}
			if ($pm == -1)
				return;

			for ($i = 0; $i < 100; $i++) {
				$dxy = $this->calcDeltaXY($pm);
				$this->xydata[$pm]['x'] = $this->xydata[$pm]['x']+$dxy[0];
				$this->xydata[$pm]['y'] = $this->xydata[$pm]['y']+$dxy[1];

				$deltam = $this->calcDeltaM($pm);
				if ($deltam < $this->EPSILON)
					break;
			}

			if ($this->adjustForGravity)
				$this->adjustForGravity();

			if ($this->exchangeVertices && $maxDeltaM < $this->EPSILON) {
				$energy = $this->calcEnergy();
				for ($i = 0; $i < $n - 1; $i++) {
					for ($j = $i + 1; $j < $n; $j++) {
						$xenergy = $this->calcEnergyIfExchanged($i, $j);
						if ($energy > $xenergy) {
							$sx = $this->xydata[$i]['x'];
							$sy = $this->xydata[$i]['y'];
							$this->xydata[$i]['x'] = $this->xydata[$j]['x'];//.setLocation(xydata[$j]);
							$this->xydata[$i]['y'] = $this->xydata[$j]['y'];
							$this->xydata[$j]['x'] = $sx;
							$this->xydata[$j]['x'] = $sy;
							return;
						}
					}
				}
			}
		}
		catch(Exception $e) {
//			fireStateChanged();
		}
	}

	/**
	 * Shift all vertices so that the center of gravity is located at
	 * the center of the screen.
	 */
	public function adjustForGravity() {
		$d = $this->size;
		$height = $d['height'];
		$width  = $d['width'];
		$gx = 0;
		$gy = 0;
		$cnt = count($this->xydata);
		for ($i = 0; $i < $cnt; $i++) {
			$gx += $this->xydata[$i]['x'];
			$gy += $this->xydata[$i]['y'];
		}
		$gx /= $cnt;
		$gy /= $cnt;
		$diffx = $width / 2 - $gx;
		$diffy = $height / 2 - $gy;
		for ($i = 0; $i < $cnt; $i++) {
            $this->xydata[$i]['x'] = $this->xydata[$i]['x']+$diffx;
            $this->xydata[$i]['y'] = $this->xydata[$i]['y']+$diffy;
		}
	}
	
/**
	 * Determines a step to new position of the vertex m.
	 */
	protected function calcDeltaXY($m) {
		$dE_dxm = 0;
		$dE_dym = 0;
		$d2E_d2xm = 0;
		$d2E_dxmdym = 0;
		$d2E_dymdxm = 0;
		$d2E_d2ym = 0;

		for ($i = 0; $i < count($this->vertices); $i++) {
			if ($i != $m) {
                
                $dist = $this->dm[$m][$i];
				$l_mi = $this->L * $dist;
				$k_mi = $this->K / ($dist * $dist);
				$dx = $this->xydata[$m]['x'] - $this->xydata[$i]['x'];
				$dy = $this->xydata[$m]['y'] - $this->xydata[$i]['y'];
				$d = sqrt($dx * $dx + $dy * $dy);
				$ddd = $d * $d * $d;

				$dE_dxm += $k_mi * (1 - $l_mi / $d) * $dx;
				$dE_dym += $k_mi * (1 - $l_mi / $d) * $dy;
				$d2E_d2xm += $k_mi * (1 - $l_mi * $dy * $dy / $ddd);
				$d2E_dxmdym += $k_mi * $l_mi * $dx * $dy / $ddd;
				$d2E_d2ym += $k_mi * (1 - $l_mi * $dx * $dx / $ddd);
			}
		}
		// d2E_dymdxm equals to d2E_dxmdym.
		$d2E_dymdxm = $d2E_dxmdym;

		$denomi = $d2E_d2xm * $d2E_d2ym - $d2E_dxmdym * $d2E_dymdxm;
		$deltaX = ($d2E_dxmdym * $dE_dym - $d2E_d2ym * $dE_dxm) / $denomi;
		$deltaY = ($d2E_dymdxm * $dE_dxm - $d2E_d2xm * $dE_dym) / $denomi;
		return array($deltaX, $deltaY);
	}

	/**
	 * Calculates the gradient of energy function at the vertex m.
	 */
	protected function calcDeltaM($m) {
		$dEdxm = 0;
		$dEdym = 0;
		for ($i = 0; $i < count($this->vertices); $i++) {
			if ($i != $m) {
                $dist = $this->dm[$m][$i];
				$l_mi = $this->L * $dist;
				$k_mi = $this->K / ($dist * $dist);

				$dx = $this->xydata[$m]['x'] - $this->xydata[$i]['x'];
				$dy = $this->xydata[$m]['y'] - $this->xydata[$i]['y'];
				$d = sqrt($dx * $dx + $dy * $dy);

				$common = $k_mi * (1 - $l_mi / $d);
				$dEdxm += $common * $dx;
				$dEdym += $common * $dy;
			}
		}
		return sqrt($dEdxm * $dEdxm + $dEdym * $dEdym);
	}

	/**
	 * Calculates the energy function E.
	 */
	protected function calcEnergy() {
		$energy = 0;
		for ($i = 0; $i < count($this->vertices)- 1; $i++) {
			for ($j = $i + 1; $j < count($this->vertices); $j++) {
                $dist = $this->dm[$i][$j];
				$l_ij = $this->L * $dist;
				$k_ij = $this->K / ($dist * $dist);
				$dx = $this->xydata[$i]['x'] - $this->xydata[$j]['x'];
				$dy = $this->xydata[$i]['y'] - $this->xydata[$j]['y'];
				$d = sqrt($dx * $dx + $dy * $dy);


				$energy += $k_ij / 2 * ($dx * $dx + $dy * $dy + $l_ij * $l_ij -
									  2 * $l_ij * $d);
			}
		}
		return $energy;
	}

	/**
	 * Calculates the energy function E as if positions of the
	 * specified vertices are exchanged.
	 */
	protected function calcEnergyIfExchanged($p, $q) {
		if ($p >= $q)
			throw new Exception("p should be < q");
		$energy = 0;		// < 0
		for ($i = 0; $i < count($this->vertices) - 1; $i++) {
			for ($j = $i + 1; $j < count($this->vertices); $j++) {
				$ii = $i;
				$jj = $j;
				if ($i == $p) $ii = $q;
				if ($j == $q) $jj = $p;

                $dist = $this->dm[$i][$j];
				$l_ij = $this->L * $dist;
				$k_ij = $this->K / ($dist * $dist);
				$dx = $this->xydata[$ii]['x'] - $this->xydata[$jj]['x'];
				$dy = $this->xydata[$ii]['y'] - $this->xydata[$jj]['y'];
				$d = sqrt($dx * $dx + $dy * $dy);
				
				$energy += $k_ij / 2 * ($dx * $dx + $dy * $dy + $l_ij * $l_ij -
									  2 * $l_ij * $d);
			}
		}
		return $energy;
	}

	public function reset() {
		$this->currentIteration = 0;
	}
	
}
