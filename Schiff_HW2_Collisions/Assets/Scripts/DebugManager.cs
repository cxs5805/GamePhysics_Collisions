using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DebugManager : MonoBehaviour
{
    // flag for whether or not to perform debug behavior
    public bool debug;
    public Material material;

    // an alternate approach: a single MTV and overlap
    public Vector3 MTV;
    public float overlapMagnitude;

    public GameObject POC;

    private const float X = 6.0f;
    private const float Y = 4.4f;

	// Use this for initialization
	void Start ()
    {
        Init();
	}
	
	// Update is called once per frame
	void Update ()
    {
        if (overlapMagnitude != float.MaxValue)
            overlapMagnitude = float.MaxValue;

        if (Input.GetKeyDown(KeyCode.Space))
        {
            debug = !debug;

            if (!debug)
                POC.transform.position = new Vector3(100.0f, 100.0f);
        }

        if (Input.GetKeyDown(KeyCode.Return))
        {
            Init();
        }
    }

    public void OnRenderObject()
    {
        if (debug)
        {
            GL.PushMatrix();

            // MTV (line drawn in red)
            material.SetPass(0);
            GL.Begin(GL.LINES);
            GL.Color(new Color(1, 0, 0, 0));
            GL.Vertex(Vector3.zero);
            GL.Vertex(MTV);
            GL.End();

            GL.PopMatrix();
        }
    }

    void Init()
    {
        // get all moving objects on scene
        GameObject[] squares = GameObject.FindGameObjectsWithTag("Square");
        GameObject[] triangles = GameObject.FindGameObjectsWithTag("Triangle");
        GameObject[] hexagons = GameObject.FindGameObjectsWithTag("Hexagon");
        //GameObject[] circles = GameObject.FindGameObjectsWithTag("Circle");

        // place randomly on scene
        Place(squares);
        Place(triangles);
        Place(hexagons);
    }

    void Place(GameObject[] objects)
    {
        for (int i = 0; i < objects.Length; i++)
        {
            // first, get movement component
            BaseMovement movement = objects[i].GetComponent<BaseMovement>();

            // random position
            objects[i].transform.position = new Vector3(Random.Range(-X, X), Random.Range(-Y, Y));

            // random velocity
            movement.velocity = new Vector3(Random.Range(-4.0f, 4.0f), Random.Range(-4.0f, 4.0f));

            // random omega
            movement.omega = Random.Range(Mathf.PI, Mathf.PI);
        }
    }
}
