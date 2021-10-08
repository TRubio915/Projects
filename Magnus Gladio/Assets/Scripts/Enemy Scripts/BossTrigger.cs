using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BossTrigger : MonoBehaviour
{


    public GameObject bossPrefab;

    public GameObject bossGate;

    public GameObject bossHealthBar;

    private GameObject bossClone;

    private float camZoomSize;

    private Camera cam;

    private void Start()
    {
        cam = Camera.main;
        camZoomSize = GameObject.FindObjectOfType<RangedAttackWithTeleportation>().maxCamZoomSize;
    }
    private void OnTriggerEnter2D(Collider2D other)
    {
        if (other.tag == "Player" && other.GetType() == typeof(BoxCollider2D) && bossClone == null)
        {

            StartCoroutine(ZoomOutCamera());
            bossClone = Instantiate(bossPrefab, bossPrefab.transform.position, bossPrefab.transform.rotation);
            
            bossHealthBar.SetActive(true);

            bossGate.GetComponent<Animator>().enabled = true;


            GetComponent<BoxCollider2D>().enabled = false;
        }

    }
    private IEnumerator ZoomOutCamera()
    {
        while (cam.orthographicSize < camZoomSize)
        {
            yield return null;

            cam.orthographicSize = Mathf.MoveTowards(cam.orthographicSize, camZoomSize, 5f * Time.unscaledDeltaTime);
        }

        GameObject.FindObjectOfType<RangedAttackWithTeleportation>().minCamZoomSize = camZoomSize;
    }
}
