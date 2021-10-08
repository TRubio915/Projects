using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Teleporter : MonoBehaviour
{
    float fade = 0;
    Material mat;

    private void Start()
    {
        mat = GetComponent<SpriteRenderer>().material;
        mat.SetFloat("_Fade", fade);
    }

    public void EnableTeleporter()
    {

      



        StartCoroutine(EnableFade());

        GetComponent<BoxCollider2D>().enabled = true;

    }

    IEnumerator EnableFade()
    {
        
        while (fade < 1)
        {
            yield return null;
            mat.SetFloat("_Fade", fade);
            fade += Time.deltaTime;

        }
    }
}