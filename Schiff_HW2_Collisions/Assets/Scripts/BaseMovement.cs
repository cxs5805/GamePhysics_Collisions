using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BaseMovement : MonoBehaviour
{
    // kinematic properties
    // transform.position
    public Vector3 velocity;
    public Vector3 acceleration;
    private const float MAX_V = 7.0f;

    // dynamic properties
    public float mass;
    public float radius;

    // angular properties
    public float theta;
    public float omega;
    public float alpha;
    private float momentOfInertia;
    private bool calculated = false;
    private const float MAX_OMEGA = 2 * Mathf.PI;

    // Camera (needed for bounds checking)
    private GameObject cameraObject;
    private Camera camera;

    // dynamic point attributes
    private PolygonCollider2D polygonCollider2D;
    private Vector3[] localPoints;
    private Vector3[] points;
    private Vector3[] localAxes;
    //private Vector3 minimumTranslationVector;
    //private float overlapMagnitude;
    private Vector3 perpendicularAxis;

    // hard-coded constant values

    // hard-coded value for floating-point epsilon
    // this is because it isn't defined in .NET Framework
    const float EPSILON = 1.19209289e-07f;

    // account for rounding errors in floating point arithmetic
    const float MARGIN_OF_ERROR = 0.01f;

    // reference to debug manager for debug
    public DebugManager debugManager;

    // Use this for initialization
    void Start()
    {
        // make sure mass and radius are both nonzero
        if (mass <= 0) mass = 1.0f;
        if (radius <= 0) radius = 0.5f;

        // define mass and radius scale
        // (use scale.x arbitrarily because we assume regular polygons)
        mass = transform.lossyScale.x;
        radius = 0.5f * transform.lossyScale.x;

        // get the camera
        cameraObject = GameObject.FindGameObjectWithTag("MainCamera");
        camera = cameraObject.GetComponent<Camera>();

        // get hard-coded local values for the collider's points' positions
        polygonCollider2D = GetComponent<PolygonCollider2D>();
        localPoints = new Vector3[polygonCollider2D.points.Length];
        for (int i = 0; i < polygonCollider2D.points.Length; i++)
        {
            localPoints[i] = new Vector3(polygonCollider2D.points[i].x, polygonCollider2D.points[i].y);
            //points[i] = new Vector3(polygonCollider2D.points[i].x, polygonCollider2D.points[i].y);
        }

        // get global equivalent
        points = new Vector3[polygonCollider2D.points.Length];
        UpdatePoints();
        //Debug.Log(UpdatePoints());

        // now that collider points are in global space,
        // we can use them to calculate moment of inertia
        CalculateMomentOfInertia();

        //Debug.Log("MTV: (" + minimumTranslationVector.x + ", " + minimumTranslationVector.y);

        // set overlap magnitude to something arbitrarily large
        //overlapMagnitude = float.MaxValue;

        // make sure debugManager is defined
        if (debugManager == null)
        {
            GameObject gameManager = GameObject.FindGameObjectWithTag("GameManager");
            debugManager = gameManager.GetComponent<DebugManager>();
        }
    }

    // Update is called once per frame
    public void Update ()
    {
        //if (minimumTranslationVector != Vector3.zero)
        //    Debug.Log("MTV: (" + minimumTranslationVector.x + ", " + minimumTranslationVector.y);

        // if need be, reset the MTV to zero and the overlap magnitude to max value
        //if (debugManager.debug)
        //{
        //minimumTranslationVector = Vector3.zero;
        //overlapMagnitude = float.MaxValue;
        //}

        // update kinematic and angular attributes
        Euler();
        CheckBounds();

        // update point attributes so collision can be done properly
        string debugStr = UpdatePoints();
        localAxes = GetLocalAxes();

        // only check collisions with like objects
        GameObject[] objectsWithSameTag = GameObject.FindGameObjectsWithTag(tag);
        foreach (GameObject tagObject in objectsWithSameTag)
        {
            BaseMovement otherMovement = tagObject.GetComponent<BaseMovement>();
            //Debug.Log(otherMovement.name);

            // ignore if the other object is actually us
            if (name == otherMovement.name)
                continue;

            //*
            // detect collision
            bool colliding = Colliding(otherMovement);

            // if there was a collision...
            if (colliding)
            {
                //Debug.Log(name + " is Colliding with " + otherMovement.name);

                // get the point of collision
                Vector3 pointOfCollision = GetPointOfCollision(otherMovement);

                GameObject POC = GameObject.FindGameObjectWithTag("POC");
                if (debugManager.debug)
                    POC.transform.position = pointOfCollision;

                // and use that to respond a la Hecker article
                float j = Respond(otherMovement, pointOfCollision);
                //Debug.Log(j);
            }
            //*/
        }

        /*
        // DEBUG
        if (Input.GetKeyDown(KeyCode.Space))
        {
            //Debug.Log(debugStr);

            //Debug.Log("local x = (" + localAxes[0].x + ", " + localAxes[0].y + ")");
            //Debug.Log("local y = (" + localAxes[1].x + ", " + localAxes[1].y + ")");
        }
        //*/
    }

    public void OnRenderObject()
    {
        /*
        if(debugManager.debug)
        {
            GL.PushMatrix();

            // MTV (line drawn in red)
            debugManager.material.SetPass(0);
            GL.Begin(GL.LINES);
            GL.Color(new Color(1, 0, 0, 0));
            GL.Vertex(Vector3.zero);
            GL.Vertex(minimumTranslationVector);
            GL.End();

            GL.PopMatrix();
        }
        //*/
    }

    void Euler()
    {
        // Euler linear
        velocity += acceleration * Time.deltaTime;

        // cap velocity
        if (velocity.magnitude >= MAX_V)
        {
            velocity.Normalize();
            velocity *= MAX_V;
        }

        transform.position += velocity * Time.deltaTime;

        // Euler angular
        omega += alpha * Time.deltaTime;

        // cap omega
        if (omega >= MAX_OMEGA)
            omega = MAX_OMEGA;

        theta += omega * Time.deltaTime;
        theta %= (Mathf.PI * 2);
        transform.rotation = Quaternion.Euler(0.0f, 0.0f, theta * Mathf.Rad2Deg);
        //Debug.Log(theta % (Mathf.PI * 2));
    }

    void CheckBounds()
    {
        // check bounds of screen
        Vector3 R = new Vector3(radius, radius);
        Vector3 minusR = camera.WorldToViewportPoint(transform.position - R);
        Vector3 plusR = camera.WorldToViewportPoint(transform.position + R);

        // manually shift positions to ensure that objects never get stuck at edges
        if (minusR.x < 0)
        {
            transform.position = new Vector3(transform.position.x + /*R.x +*/ 0.001f, transform.position.y);
        }
        else if (plusR.x > 1)
        {
            transform.position = new Vector3(transform.position.x - /*R.x -*/ 0.001f, transform.position.y);
        }

        if (minusR.y < 0)
        {
            transform.position = new Vector3(transform.position.x, transform.position.y + /*R.y +*/ 0.001f);
        }
        else if (plusR.y > 1)
        {
            transform.position = new Vector3(transform.position.x, transform.position.y - /*R.y +*/ 0.001f);
        }

        if (minusR.x < 0 || plusR.x > 1)
        {
            velocity.x = -velocity.x;
            acceleration.x = -acceleration.x;
        }

        if (minusR.y < 0 || plusR.y > 1)
        {
            velocity.y = -velocity.y;
            acceleration.y = -acceleration.y;
        }
    }

    void CalculateMomentOfInertia()
    {
        if (!calculated)
        {
            // moment of inertia only needs to be calculated once
            calculated = true;
            float L = 0.0f;
            switch (tag)
            {
                case "Square":
                    //Debug.Log("Square");

                    // for a square, the radius member data represents half of L
                    L = 2 * radius;

                    // thus, I_square =  1/12 M L^2
                    // assuming uniform density throughout square
                    momentOfInertia = 0.08333f * mass * L * L; 
                    break;
                case "Triangle":
                    //Debug.Log("Triangle");

                    // get side length L
                    L = Vector3.Distance(points[1], points[2]);

                    // I_triangle = 1/12 M L^2
                    momentOfInertia = 0.08333f * mass * L * L;
                    break;
                case "Circle":
                    //Debug.Log("Circle");

                    // I_circle = 1/2 M R^2
                    momentOfInertia = 0.5f * mass * radius * radius;
                    break;
                case "Hexagon":
                    //Debug.Log("Hexagon");

                    // regular hexagon can be composed of 6 equilateral triangles
                    // I_hexagon = 5/6 M L^2 (?)
                    L = Vector3.Distance(points[0], points[1]);
                    momentOfInertia = 0.83333f * mass * L * L;
                    break;
                default:
                    break;
            }
        }
    }

    string UpdatePoints()
    {
        int numPoints = localPoints.Length;
        string debugStr = "global center: ( " + transform.position.x + ", " + transform.position.y + ")\n";
        for (int i = 0; i < numPoints; i++)
        {
            // local
            //debugStr += "local point " + i + ": (" + localPoints[i].x + ", " + localPoints[i].y + ")\n";

            // get global points
            //Vector3 localPoint = new Vector3(localPoints[i].x, localPoints[i].y);
            points[i] = transform.TransformPoint(localPoints[i]);
            
            debugStr += "global point " + i + ": (" + points[i].x + ", " + points[i].y + ")\n";
        }
        return debugStr;
    }

    Vector3[] GetLocalAxes()
    {
        Vector3[] axes = new Vector3[2];

        Vector3 startX = Vector3.zero;
        Vector3 endX = Vector3.zero;
        Vector3 startY = Vector3.zero;
        Vector3 endY = Vector3.zero;

        switch(tag)
        {
            case "Square":
                //Debug.Log("Square");

                startX = points[0];
                endX = points[3];
                startY = startX;
                endY = points[1];

                break;
            case "Triangle":
                startX = points[1];
                endX = points[2];
                startY = new Vector3(points[0].x, points[1].y);
                endY = points[0];
                break;
            case "Hexagon":
                startX = points[1];
                endX = points[5];
                startY = startX;
                endX = points[2];
                break;
            default: break;
        }

        // index 0 is x
        axes[0] = GetLocalAxis(startX, endX);

        // index 1 is y
        axes[1] = GetLocalAxis(startY, endY);

        return axes;
    }

    Vector3 GetLocalAxis(Vector3 start, Vector3 end)
    {
        return Vector3.Normalize(end - start);
    }

    bool Colliding(BaseMovement other)
    {
        // SAT
        // get my local X and Y
        Vector3 myX = localAxes[0];
        Vector3 myY = localAxes[1];

        // get other object's local X and Y
        Vector3 otherX = other.localAxes[0];
        Vector3 otherY = other.localAxes[1];

        // now we can actually check collisions along the 4 axes
        // my x
        if (!CollidingOnAxis(myX, other))
        {
            return false;
        }

        // my y
        if (!CollidingOnAxis(myY, other))
        {
            return false;
        }

        // other x
        if (!CollidingOnAxis(otherX, other))
        {
            return false;
        }

        // other y
        if (!CollidingOnAxis(otherY, other))
        {
            return false;
        }

        // all axes failed => they were colliding
        // at this point, the MTV is up to date
        // so make sure that MTV always faces us, the current object
        Vector3 otherToMe = transform.position - other.transform.position;
        //if (Vector3.Dot(minimumTranslationVector, otherToMe) < 0.0f)
        //    minimumTranslationVector *= -1.0f;
        if (Vector3.Dot(debugManager.MTV, otherToMe) < 0.0f)
            debugManager.MTV *= -1.0f;

        return true;
    }

    bool CollidingOnAxis(Vector3 axis, BaseMovement other)
    {
        // we need to project my global points on the axis
        float[] myPointsOnAxis = new float[points.Length];
        for (int i = 0; i < points.Length; i++)
        {
            myPointsOnAxis[i] = Vector3.Dot(points[i], axis);
        }

        // do the same for the other's global points
        float[] otherPointsOnAxis = new float[other.points.Length];
        for (int i = 0; i < other.points.Length; i++)
        {
            otherPointsOnAxis[i] = Vector3.Dot(other.points[i], axis);
        }

        // find my max and min values along axis
        float myMaxOnAxis = float.MinValue;
        float myMinOnAxis = float.MaxValue;
        for (int i = 0; i < points.Length; i++)
        {
            if (myPointsOnAxis[i] > myMaxOnAxis)
                myMaxOnAxis = myPointsOnAxis[i];
            else if (myPointsOnAxis[i] < myMinOnAxis)
                myMinOnAxis = myPointsOnAxis[i];
        }

        // now do the same for the other
        float otherMaxOnAxis = float.MinValue;
        float otherMinOnAxis = float.MaxValue;
        for (int i = 0; i < other.points.Length; i++)
        {
            if (otherPointsOnAxis[i] > otherMaxOnAxis)
                otherMaxOnAxis = otherPointsOnAxis[i];
            else if (otherPointsOnAxis[i] < otherMinOnAxis)
                otherMinOnAxis = otherPointsOnAxis[i];
        }

        // time to test the collision on the axis itself
        if (myMaxOnAxis < otherMinOnAxis || otherMaxOnAxis < myMinOnAxis)
        {
            return false;
        }
        else
        {
            // if they are colliding, update the minimum translation vector
            //Debug.Log("There is a collision on this axis");

            // get the 2 overlaps about the axis
            float outerOverlap = otherMaxOnAxis - myMinOnAxis;
            float innerOverlap = myMaxOnAxis - otherMinOnAxis;

            // get the smaller of the 2
            float minOverlap = outerOverlap;
            if (innerOverlap < outerOverlap)
                minOverlap = innerOverlap;

            // if the smallest overlap in this collision is the smallest
            // overlap we have seen thus far,
            //if (minOverlap < overlapMagnitude)
            if (minOverlap < debugManager.overlapMagnitude)
            {
                //overlapMagnitude = minOverlap;
                //minimumTranslationVector = axis;

                debugManager.overlapMagnitude = minOverlap;
                debugManager.MTV = axis;

                //Debug.Log("MTV set! (" + minimumTranslationVector.x + ", " + minimumTranslationVector.y + ")");
            }
            
        }
        return true;
    }

    Vector3 GetPointOfCollision(BaseMovement other)
    {
        // based on ATLAS Physics Level 8: Determining Point of Collision 2D (Convex Hull)
        // code: https://github.com/IGME-RIT/physics-determiningCollisionPoint-ConvexHull-2D

        List<Vector3> myClosestPoints = new List<Vector3>();
        List<Vector3> otherClosestPoints = new List<Vector3>();

        // get the point that is the least in the direction of my MTV
        // (this helper method will also update currentPoints properly)
        float minAlongMyMTV = MinAlongMyMTV(myClosestPoints);

        // if there's only one closest point to us, then just return that
        if (myClosestPoints.Count == 1)
            return myClosestPoints[0];

        // get the point that is the greatest in the direction of the other MTV
        float maxAlongOtherMTV = MaxAlongOtherMTV(other, otherClosestPoints);

        // again, return sole element if there's only one closest point to other
        if (otherClosestPoints.Count == 1)
            return otherClosestPoints[0];

        // since there are multiple points, we must check which one
        // is closest on an axis perpendicular to our MTV
        //perpendicularAxis = new Vector3(-minimumTranslationVector.y, minimumTranslationVector.x);
        perpendicularAxis = new Vector3(-debugManager.MTV.y, debugManager.MTV.x);
        perpendicularAxis.Normalize();

        // combine the two lists of points into one
        List<Vector3> closestPoints = new List<Vector3>();

        // add from my points
        for (int i = 0; i < myClosestPoints.Count; i++)
            closestPoints.Add(myClosestPoints[i]);
        // and add from other points
        for (int i = 0; i < otherClosestPoints.Count; i++)
            closestPoints.Add(otherClosestPoints[i]);

        // determine which points are the min and max along edge axis
        float minTemp = Vector3.Dot(closestPoints[0], perpendicularAxis);
        float maxTemp = minTemp;

        int minIndex = 0;
        int maxIndex = 0;

        for (int i = 0; i < closestPoints.Count; i++)
        {
            // project point onto edge axis
            float dot = Vector3.Dot(closestPoints[i], perpendicularAxis);

            // update min and max
            if (dot < minTemp)
            {
                minTemp = dot;
                minIndex = i;
            }
            if (dot > maxTemp)
            {
                maxTemp = dot;
                maxIndex = i;
            }
        }

        // remove max and min from the combined list
        closestPoints.RemoveAt(minIndex);
        if (minIndex < maxIndex)
            maxIndex--;
        closestPoints.RemoveAt(maxIndex);

        // average the first two points after max and min are gone
        // the result should be close enough to what the closest point would be!
        Vector3 theClosestPoint = (closestPoints[0] + closestPoints[1]) * 0.5f;
        return theClosestPoint;
    }

    float MinAlongMyMTV(List<Vector3> closestPoints)
    {
        // initially, assume first point is our minimum along MTV
        // i.e. since the MTV faces the current object,
        // the minimum is the closest point along the MTV to the OTHER object
        Vector3 currentPoint = points[0];
        //float currentMin = Vector3.Dot(currentPoint, minimumTranslationVector);
        float currentMin = Vector3.Dot(currentPoint, debugManager.MTV);

        // keep track of all of the points tied for the min
        //closestPoints = new List<Vector3>(); 

        // loop thru the object's points (already in global space)
        int numPoints = points.Length;
        for (int i = 0; i < numPoints; i++)
        {
            // project current point onto MTV
            // to get the distance of the point towards the other object
            currentPoint = points[i];
            //float dot = Vector3.Dot(currentPoint, minimumTranslationVector);
            float dot = Vector3.Dot(currentPoint, debugManager.MTV);

            // if the current point is same distance away as already known min
            if (Mathf.Abs(dot - currentMin) < EPSILON + MARGIN_OF_ERROR)
            {
                // add the new point to our list of closest points
                closestPoints.Add(currentPoint);
            }
            else if (dot < currentMin - EPSILON)
            {
                // new min
                currentMin = dot;

                // clear the list because it's the first point this far back on MTV
                // so there are no ties at the moment
                closestPoints.Clear();
                closestPoints.Add(currentPoint);
            }

        }

        return currentMin;
    }

    float MaxAlongOtherMTV(BaseMovement other, List<Vector3> closestPoints)
    {
        // this function is almost identitical to MinAlongMyMTV
        // with two key differences:
        // 1. we track a MAX instead of a MIN
        // 2. it's along the MTV of the OTHER object
        Vector3 currentPoint = other.points[0];
        //float currentMax = Vector3.Dot(currentPoint, other.minimumTranslationVector);
        float currentMax = Vector3.Dot(currentPoint, debugManager.MTV);

        // variable to help debug
        //Vector3 mtv = other.minimumTranslationVector;
        Vector3 mtv = debugManager.MTV;
        //mtv = minimumTranslationVector;

        int numPoints = other.points.Length;
        for (int i = 0; i < numPoints; i++)
        {
            currentPoint = other.points[i];

            //float dot = Vector3.Dot(currentPoint, other.minimumTranslationVector);
            float dot = Vector3.Dot(currentPoint, mtv);

            if (Mathf.Abs(dot - currentMax) < EPSILON + MARGIN_OF_ERROR)
            {
                closestPoints.Add(currentPoint);
            }
            else if (dot > currentMax + EPSILON)
            {
                currentMax = dot;
                closestPoints.Clear();
                closestPoints.Add(currentPoint);
            }
        }

        return currentMax;
    }

    float Respond(BaseMovement other, Vector3 pointOfCollision)
    {
        float j = 0.0f;

        // get relative velocity
        // calculate vectors from both centers of mass to the point of collision
        Vector3 rAP = pointOfCollision - transform.position;
        Vector3 rBP = pointOfCollision - other.transform.position;

        // calculate relative velocities of the point of collision on both objects
        Vector3 vAP = omega * rAP;
        Vector3 vBP = other.omega * rBP;

        // get the total relative velocity
        Vector3 relativeVelocity = vAP - vBP;
        //Debug.Log("rel v = (" + relativeVelocity.x + ", " + relativeVelocity.y);

        // numerator

        // energy is conserved
        float e = 1.0f;
        float eFactor = -1 * (1.0f + e);
        Vector3 vWithFactor = eFactor * relativeVelocity;
        //float numerator = Vector3.Dot(vWithFactor, minimumTranslationVector);
        float numerator = Vector3.Dot(vWithFactor, debugManager.MTV);
        //Debug.Log("numerator = " + numerator);

        // denominator
        // linear component
        float massFactor = (1 / mass) + (1 / other.mass);
        //Vector3 nWithFactor = minimumTranslationVector * massFactor;
        Vector3 nWithFactor = debugManager.MTV * massFactor;
        //float linearDenominator = Vector3.Dot(minimumTranslationVector, nWithFactor);
        float linearDenominator = Vector3.Dot(debugManager.MTV, nWithFactor);
        //Debug.Log("lin = " + linearDenominator);

        // angular component
        //float myDot = Vector3.Dot(perpendicularAxis, minimumTranslationVector);
        //float myDot = Vector3.Dot(rAP, minimumTranslationVector);
        float myDot = Vector3.Dot(rAP, debugManager.MTV);
        float myAngular = (myDot * myDot) / momentOfInertia;

        //float otherDot = Vector3.Dot(other.perpendicularAxis, other.minimumTranslationVector);
        //float otherDot = Vector3.Dot(rBP, other.minimumTranslationVector);
        float otherDot = Vector3.Dot(rBP, debugManager.MTV);
        float otherAngular = (otherDot * otherDot) / momentOfInertia;

        float angularDenominator = myAngular + otherAngular;
        //Debug.Log("ang = " + angularDenominator);

        // linear + angular
        float denominator = linearDenominator + angularDenominator;
        //Debug.Log("denom = " + denominator);

        // j = numerator / denominator
        j = numerator / denominator;
        //Debug.Log("numerator / denominator = " + numerator + "/" + denominator);
        //Debug.Log("j = " + j);

        // now, the collision response itself
        // we only do this for the current object so as to not apply the correct impulse one too many times
        // apply linear impulse
        float jPerM = j / mass;
        //velocity = velocity + (jPerM * minimumTranslationVector);
        velocity = velocity + (jPerM * debugManager.MTV);

        // apply angular impulse
        //Vector3 jN = j * minimumTranslationVector;
        Vector3 jN = j * debugManager.MTV;
        float perpDot = Vector3.Dot(perpendicularAxis, jN);
        omega = omega + (perpDot / momentOfInertia);

        // return j value for debugging
        return j;
    }
}
